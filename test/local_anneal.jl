# Local Annealing ZIT — adaptive price beliefs (§17) + signal-gated annealing (§13 style).
# Fully decentralized: no p*, D, S, or q* ever observed.
#
# Agents observe ONLY: own WTP/WTAC array + own realized k_i per period.
#
# State per agent: (p̂_i, λ_i)
#
# ── Quantity draw (mixed strategy) ───────────────────────────────────────────
#   With prob λ_i:   BiGeo over active units {WTP ≥ p̂_i}, p̂_i as lower bound  [explore]
#   With prob 1-λ_i: Standard ZIT — draw uniformly from [0, n_total]            [exploit]
#
# ── λ_i update (marginal surplus relative to p̂_i) ────────────────────────────
#   WTP[k] ≥ p̂_i (positive surplus at current belief):  λ *= (1−β)    [cool → exploit]
#   WTP[k] < p̂_i (negative surplus — standard ZIT over-shot): λ += (1−λ)*α [heat → explore]
#   k = 0:                                                hold
#
#   Why this equilibrates: BiGeo always draws WTP[k] ≥ p̂_i (positive → cool).
#   Standard ZIT sometimes draws beyond n_active(p̂_i) → WTP[k] < p̂_i (negative → heat).
#   Cooling / heating balance at λ_eq ≈ φ⁻² ≈ 0.38 (same golden-ratio result as §13,
#   now with p̂_i replacing p* in the surplus comparison).
#
# ── p̂_i update (local price discovery, EMA rate η) ──────────────────────────
#   Only updates from BiGeo draws — standard ZIT rationing is too noisy:
#   when all agents over-demand randomly, D>>S yields WTP[1] (highest unit) as
#   the signal, pushing p̂ far above p*.
#   BiGeo rationed (k < d):   p̂ += η*(WTP[k]    − p̂)  [raise toward marginal WTP received]
#   BiGeo at cap (d = na):    p̂ += η*(WTP[na+1] − p̂)  [lower toward first excluded unit]
#   Otherwise:                hold
#
# ── Why this works ────────────────────────────────────────────────────────────
#   §13 needed p* for BiGeo lower bound and surplus comparison.
#   §17 gave local p̂ but had equilibrium drift as "at cap" fires every turn at WGE.
#   Here, p̂ replaces p* in both roles. λ equilibrates at ~0.38 (not 0) because
#   standard ZIT's random over-demand provides the heat signal that keeps λ > 0,
#   maintaining stochastic pressure that prevents false-equilibrium lock-in (§16).
#
# Init: p̂_i = mean(own WTP), p̂_j = mean(own WTAC), λ_i = λ_j = 1.0 (pure BiGeo)

using Random
using Statistics
using Printf
using DiscreteMarket

const ETA   = 0.20   # price belief EMA rate
const ALPHA = 0.05   # λ heat rate when rationed
const BETA  = 0.05   # λ cool rate when at cap
const T_TURNS = 300

# ── Market-level surplus ──────────────────────────────────────────────────────

function surplus_f64(m::GoodMarket, q::Int) :: Float64
    q <= 0 && return 0.0
    all_wtp  = sort([Float64(v) for d in m.consumers for v in d.wtp], rev=true)
    all_wtac = sort([Float64(c) for f in m.firms     for c in f.wtac])
    n = min(q, length(all_wtp), length(all_wtac))
    n == 0 && return 0.0
    sum(all_wtp[1:n]) - sum(all_wtac[1:n])
end

# ── Pro-rata realized trade ───────────────────────────────────────────────────

@inline realized_c(draw, D, S) = D == 0 ? 0 : (D <= S ? draw : floor(Int, draw * S / D))
@inline realized_f(draw, D, S) = S == 0 ? 0 : (S <= D ? draw : floor(Int, draw * D / S))

# ── Per-agent data ────────────────────────────────────────────────────────────

struct LocalConsumer
    wtp :: Vector{Float64}   # all valuations, sorted descending
end

struct LocalFirm
    wtac :: Vector{Float64}  # all costs, sorted ascending
end

LocalConsumer(d::ConsumerDemand) = LocalConsumer(sort(Float64.(d.wtp), rev=true))
LocalFirm(f::FirmSupply)         = LocalFirm(sort(Float64.(f.wtac)))

@inline n_active_c(c::LocalConsumer, p::Float64) = count(v -> v >= p, c.wtp)
@inline n_active_f(f::LocalFirm,     p::Float64) = count(v -> v <= p, f.wtac)

# ── BiGeo draw over active units using p̂ as lower bound ──────────────────────

@inline function draw_bigeo_c(c::LocalConsumer, p_hat::Float64, na::Int,
                               rng::AbstractRNG) :: Int
    na == 0 && return 0
    log_lo = log(max(p_hat, 1e-10))
    log_hi = log(c.wtp[1])                 # max WTP (= highest active since wtp sorted desc)
    log_hi <= log_lo + 1e-12 && return na  # p̂ ≈ max WTP → buy all active
    pv = exp(log_lo + rand(rng) * (log_hi - log_lo))
    # count(WTP ≥ pv) in sorted-descending array
    searchsortedlast(c.wtp, pv, rev=true)
end

@inline function draw_bigeo_f(f::LocalFirm, p_hat::Float64, na::Int,
                               rng::AbstractRNG) :: Int
    na == 0 && return 0
    log_lo = log(f.wtac[1])                # min WTAC (cheapest, always active if na>0)
    log_hi = log(max(p_hat, f.wtac[1] + 1e-10))
    log_hi <= log_lo + 1e-12 && return na  # p̂ ≈ min WTAC → supply 1 unit
    pv = exp(log_lo + rand(rng) * (log_hi - log_lo))
    # count(WTAC ≤ pv) in sorted-ascending array
    searchsortedlast(f.wtac, pv)
end

# ── Mixed draw: BiGeo with prob λ, Standard ZIT with prob 1-λ ────────────────
# Returns (d, na, from_bigeo): demand submitted, n_active at current p̂, draw type.
# Standard ZIT can demand beyond n_active — this is what generates the heat signal for λ.

@inline function draw_c(c::LocalConsumer, p_hat::Float64, λ::Float64,
                         rng::AbstractRNG) :: Tuple{Int,Int,Bool}
    na = n_active_c(c, p_hat)
    rand(rng) < λ ? (draw_bigeo_c(c, p_hat, na, rng), na, true) :
                    (rand(rng, 0:length(c.wtp)),        na, false)
end

@inline function draw_f(f::LocalFirm, p_hat::Float64, λ::Float64,
                         rng::AbstractRNG) :: Tuple{Int,Int,Bool}
    na = n_active_f(f, p_hat)
    rand(rng) < λ ? (draw_bigeo_f(f, p_hat, na, rng), na, true) :
                    (rand(rng, 0:length(f.wtac)),       na, false)
end

# ── Price belief update ───────────────────────────────────────────────────────
# "At cap, lower p̂" only fires for BiGeo draws: standard ZIT hitting exactly na is noise.

@inline function update_p_c(p::Float64, wtp::Vector{Float64},
                              d::Int, k::Int, na::Int, from_bigeo::Bool) :: Float64
    !from_bigeo && return p  # ZIT draws are too noisy for p̂: D>>S makes WTP[k] >> p*
    if k < d
        pt = k > 0 ? wtp[k] : wtp[1]
        return (1.0 - ETA) * p + ETA * pt
    elseif d == na && na < length(wtp)
        return (1.0 - ETA) * p + ETA * wtp[na + 1]
    end
    return p
end

@inline function update_p_f(p::Float64, wtac::Vector{Float64},
                              d::Int, k::Int, na::Int, from_bigeo::Bool) :: Float64
    !from_bigeo && return p  # ZIT draws are too noisy for p̂: S>>D makes WTAC[k] << p*
    if k < d
        pt = k > 0 ? wtac[k] : wtac[1]
        return (1.0 - ETA) * p + ETA * pt
    elseif d == na && na < length(wtac)
        return (1.0 - ETA) * p + ETA * wtac[na + 1]
    end
    return p
end

# ── λ update (marginal surplus relative to p̂) ────────────────────────────────
# BiGeo draws always yield WTP[k] ≥ p̂ (positive surplus) → cool.
# Standard ZIT draws beyond n_active yield WTP[k] < p̂ (negative surplus) → heat.
# The heat/cool balance stabilises λ at λ_eq ≈ φ⁻² ≈ 0.38 (golden-ratio, as in §13).

@inline function update_λ_c(λ::Float64, wtp::Vector{Float64},
                              k::Int, p_hat::Float64) :: Float64
    k == 0 && return λ
    wtp[k] >= p_hat ? λ * (1.0 - BETA) : λ + (1.0 - λ) * ALPHA
end

@inline function update_λ_f(λ::Float64, wtac::Vector{Float64},
                              k::Int, p_hat::Float64) :: Float64
    k == 0 && return λ
    wtac[k] <= p_hat ? λ * (1.0 - BETA) : λ + (1.0 - λ) * ALPHA
end

# ── Simulate ──────────────────────────────────────────────────────────────────

function simulate_local_anneal(consumers::Vector{LocalConsumer},
                                firms::Vector{LocalFirm},
                                m::GoodMarket, s_wge::Float64, q_star::Int;
                                seed::Int=0) :: Tuple{Vector{Float64},Vector{Float64},
                                                      Vector{Float64},Vector{Float64}}

    rng = MersenneTwister(seed)
    nc, nf = length(consumers), length(firms)

    p_hat_c  = [mean(c.wtp)  for c in consumers]
    p_hat_f  = [mean(f.wtac) for f in firms]
    lambda_c = ones(Float64, nc)
    lambda_f = ones(Float64, nf)

    effs   = Vector{Float64}(undef, T_TURNS)
    q_frac = Vector{Float64}(undef, T_TURNS)
    λ_c_t  = Vector{Float64}(undef, T_TURNS)
    λ_f_t  = Vector{Float64}(undef, T_TURNS)

    for t in 1:T_TURNS
        # Each draw returns (d, na, from_bigeo)
        drc = [draw_c(consumers[i], p_hat_c[i], lambda_c[i], rng) for i in 1:nc]
        drf = [draw_f(firms[j],     p_hat_f[j], lambda_f[j], rng) for j in 1:nf]

        D = sum(r[1] for r in drc; init=0)
        S = sum(r[1] for r in drf; init=0)
        q = min(D, S)

        effs[t]   = surplus_f64(m, q) / s_wge
        q_frac[t] = q_star > 0 ? q / q_star : 1.0
        λ_c_t[t]  = mean(lambda_c)
        λ_f_t[t]  = mean(lambda_f)

        ks_c = [realized_c(drc[i][1], D, S) for i in 1:nc]
        ks_f = [realized_f(drf[j][1], D, S) for j in 1:nf]

        for i in 1:nc
            d, na, bg = drc[i]; k = ks_c[i]
            old_p = p_hat_c[i]
            p_hat_c[i]  = update_p_c(old_p, consumers[i].wtp, d, k, na, bg)
            lambda_c[i] = update_λ_c(lambda_c[i], consumers[i].wtp, k, old_p)
        end
        for j in 1:nf
            d, na, bg = drf[j]; k = ks_f[j]
            old_p = p_hat_f[j]
            p_hat_f[j]  = update_p_f(old_p, firms[j].wtac, d, k, na, bg)
            lambda_f[j] = update_λ_f(lambda_f[j], firms[j].wtac, k, old_p)
        end
    end

    effs, q_frac, λ_c_t, λ_f_t
end

# ── Baseline: standard ZIT (no p* needed) ────────────────────────────────────

function baseline_std(consumers, firms, m, s_wge, rng)
    nc, nf = length(consumers), length(firms)
    mean(begin
        dc = sum(rand(rng, 0:length(consumers[i].wtp)) for i in 1:nc; init=0)
        ds = sum(rand(rng, 0:length(firms[j].wtac))    for j in 1:nf; init=0)
        surplus_f64(m, min(dc, ds)) / s_wge
    end for _ in 1:T_TURNS)
end

# ── Run across markets ────────────────────────────────────────────────────────

function run_local_anneal(; nc, nf, n_markets, seed_base, label,
                           max_units=5, report_turns=[1,2,3,5,10,20,50,100,200,300])

    sum_eff    = zeros(T_TURNS); sum_eff2   = zeros(T_TURNS)
    sum_qfrac  = zeros(T_TURNS)
    sum_λc     = zeros(T_TURNS); sum_λf     = zeros(T_TURNS)
    n_analyzed = 0
    std_effs   = Float64[]

    for seed in 1:n_markets
        rng = MersenneTwister(seed + seed_base)
        m, r = generate_good_market(rng;
            good=1, n_consumers=nc, n_firms=nf,
            max_units=max_units, Q=100)
        r.cleared || continue
        p_star = r.price
        q_star = aggregate_demand(m, p_star)
        s_wge  = surplus_f64(m, q_star)
        s_wge > 0 || continue

        consumers = [LocalConsumer(d) for d in m.consumers]
        firms     = [LocalFirm(f)     for f in m.firms]

        effs, qf, λc, λf = simulate_local_anneal(consumers, firms, m, s_wge, q_star;
                                                   seed=seed + seed_base + 500_000)
        for t in 1:T_TURNS
            sum_eff[t]   += effs[t];  sum_eff2[t]  += effs[t]^2
            sum_qfrac[t] += qf[t]
            sum_λc[t]    += λc[t];    sum_λf[t]    += λf[t]
        end
        n_analyzed += 1

        rng2 = MersenneTwister(seed + seed_base + 99_000)
        push!(std_effs, baseline_std(consumers, firms, m, s_wge, rng2))
    end

    n = n_analyzed
    println("=" ^ 70)
    println(label)
    @printf("  n=%d  T=%d  η=%.2f  α=%.2f  β=%.2f\n\n", n, T_TURNS, ETA, ALPHA, BETA)
    @printf("  Baseline — Static Std:   mean eff = %.4f\n", mean(std_effs))
    @printf("  Target   — WGE:          mean eff = 1.0000\n")
    println()
    println("  Turn │  q/q*    Eff_mean  Eff_std   λ_c     λ_f    vs Std   vs WGE")
    println("  ─────┼──────────────────────────────────────────────────────────────")
    for t in report_turns
        t > T_TURNS && continue
        eff = sum_eff[t] / n;  qf  = sum_qfrac[t] / n
        sd  = sqrt(max(0.0, sum_eff2[t]/n - eff^2))
        lc  = sum_λc[t] / n;   lf  = sum_λf[t] / n
        @printf("  t=%3d │  %.4f  %.4f    %.4f   %.3f  %.3f   %+.4f  %+.4f\n",
                t, qf, eff, sd, lc, lf, eff - mean(std_effs), eff - 1.0)
    end
    println()
    fe = sum_eff[T_TURNS]/n; fq = sum_qfrac[T_TURNS]/n
    lc = sum_λc[T_TURNS]/n;  lf = sum_λf[T_TURNS]/n
    @printf("  Final (t=%d): q/q*=%.4f  eff=%.4f  std=%.4f  λ_c=%.3f  λ_f=%.3f  vs WGE=%+.4f\n\n",
            T_TURNS, fq, fe, sqrt(max(0.0,sum_eff2[T_TURNS]/n - fe^2)), lc, lf, fe - 1.0)

    (sum_eff ./ n, sum_qfrac ./ n, n)
end

# ── Main ──────────────────────────────────────────────────────────────────────

println("LOCAL ANNEALING ZIT — adaptive price beliefs + marginal-surplus annealing (§18v2)")
println("  Agents observe ONLY: own WTP/WTAC + own realized k_i")
println("  State: (p̂_i, λ_i) per agent  [no prev_k needed]")
println()
println("  Quantity:  prob λ   → BiGeo over {WTP ≥ p̂_i}, p̂_i as lower bound  [explore]")
println("             prob 1-λ → Standard ZIT: uniform draw [0, n_total]         [exploit]")
println("  Price:     BiGeo rationed (k < d) → p̂ += η*(WTP[k] - p̂)        [raise; BiGeo only]")
println("             BiGeo at cap (d=na, k=d) → p̂ += η*(WTP[na+1] - p̂)  [lower; BiGeo only]")
println("             (ZIT rationing excluded: D>>S makes signal WTP[k]>>p*, too noisy)")
println("  Lambda:    WTP[k] ≥ p̂ (positive surplus) → λ *= (1-β)         [cool → exploit]")
println("             WTP[k] < p̂ (negative surplus)  → λ += (1-λ)*α       [heat → explore]")
println("             k=0                              → hold")
println("  Init:      p̂ = mean(WTP),  λ = 1.0  (pure BiGeo)")
println("  λ_eq theory: BiGeo always positive → cools; std ZIT over-shoots → heats.")
println("               Balance at λ_eq ≈ φ⁻² ≈ 0.38, same as §13 with p̂ replacing p*.")
println()

RT = [1, 2, 3, 5, 10, 20, 50, 100, 200, 300]

# Small / mixed markets
let
    SEED_BASE = 80_000
    sum_eff   = zeros(T_TURNS); sum_eff2  = zeros(T_TURNS)
    sum_qfrac = zeros(T_TURNS)
    sum_λc    = zeros(T_TURNS); sum_λf    = zeros(T_TURNS)
    n         = 0
    std_effs  = Float64[]

    for seed in 1:1_000
        rng  = MersenneTwister(seed + SEED_BASE)
        n_c  = rand(rng, 1:8); n_f = rand(rng, 1:8); max_u = rand(rng, 1:6)
        m, r = generate_good_market(rng;
            good=1, n_consumers=n_c, n_firms=n_f, max_units=max_u, Q=100)
        r.cleared || continue
        p_star = r.price
        q_star = aggregate_demand(m, p_star)
        s_wge  = surplus_f64(m, q_star)
        s_wge > 0 || continue

        consumers = [LocalConsumer(d) for d in m.consumers]
        firms     = [LocalFirm(f)     for f in m.firms]

        effs, qf, λc, λf = simulate_local_anneal(consumers, firms, m, s_wge, q_star;
                                                   seed=seed + SEED_BASE + 500_000)
        for t in 1:T_TURNS
            sum_eff[t]  += effs[t]; sum_eff2[t] += effs[t]^2
            sum_qfrac[t] += qf[t]
            sum_λc[t]   += λc[t];   sum_λf[t]  += λf[t]
        end
        n += 1

        rng2 = MersenneTwister(seed + SEED_BASE + 99_000)
        push!(std_effs, baseline_std(consumers, firms, m, s_wge, rng2))
    end

    println("=" ^ 70)
    println("Small/mixed markets  (1–8C, 1–8F, 1–6 max_units)  n=$n / 1000")
    @printf("  T=%d  η=%.2f  α=%.2f  β=%.2f\n\n", T_TURNS, ETA, ALPHA, BETA)
    @printf("  Baseline — Static Std:   mean eff = %.4f\n", mean(std_effs))
    println()
    println("  Turn │  q/q*    Eff_mean  Eff_std   λ_c     λ_f    vs Std   vs WGE")
    println("  ─────┼──────────────────────────────────────────────────────────────")
    for t in RT
        t > T_TURNS && continue
        eff = sum_eff[t]/n; qf = sum_qfrac[t]/n
        sd  = sqrt(max(0.0, sum_eff2[t]/n - eff^2))
        lc = sum_λc[t]/n;  lf = sum_λf[t]/n
        @printf("  t=%3d │  %.4f  %.4f    %.4f   %.3f  %.3f   %+.4f  %+.4f\n",
                t, qf, eff, sd, lc, lf, eff - mean(std_effs), eff - 1.0)
    end
    println()
    fe = sum_eff[T_TURNS]/n; fq = sum_qfrac[T_TURNS]/n
    lc = sum_λc[T_TURNS]/n;  lf = sum_λf[T_TURNS]/n
    @printf("  Final (t=%d): q/q*=%.4f  eff=%.4f  std=%.4f  λ_c=%.3f  λ_f=%.3f  vs WGE=%+.4f\n\n",
            T_TURNS, fq, fe, sqrt(max(0.0,sum_eff2[T_TURNS]/n - fe^2)), lc, lf, fe - 1.0)
end

run_local_anneal(nc=50, nf=50, n_markets=300, seed_base=90_000,
                 label="100 agents (50C + 50F)  n=300 markets",
                 report_turns=RT)

run_local_anneal(nc=250, nf=250, n_markets=150, seed_base=91_000,
                 label="500 agents (250C + 250F)  n=150 markets",
                 report_turns=RT)
