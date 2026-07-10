# §19 Local Elasticity Annealing — BiGeo + ratchet with elasticity-gap λ signal.
# Fully decentralized: no p*, D, S, or q* ever observed.
#
# Agents observe ONLY: own WTP/WTAC array + own realized k_i per period.
#
# State per agent: (p̂_i, λ_i, prev_k_i)
#
# ── Quantity draw (mixed strategy) ───────────────────────────────────────────
#   With prob λ_i:   BiGeo over {WTP ≥ p̂_i}, p̂_i as lower bound   [explore]
#   With prob 1-λ_i: Ratchet = min(n_active(p̂_i), prev_k_i + 1)   [exploit]
#   Both are bounded by n_active(p̂_i) so neither over-demands beyond the active set.
#
# ── p̂_i update (§17 three-case EMA, valid for both draw types) ──────────────
#   Both draw mechanisms are bounded by na(p̂_i), so rationing and at-cap signals
#   are unambiguous (unlike ZIT which could over-demand past na):
#   Rationed (k < d):       p̂ += η*(WTP[k]    − p̂)   [raise]
#   At cap (d = na, k = d): p̂ += η*(WTP[na+1] − p̂)   [lower]
#   Otherwise:              hold
#
# ── λ_i update (local elasticity gap signal) ─────────────────────────────────
#   gap  = WTP[k_i] − WTP[k_i+1]              (forward drop in own WTP at traded unit)
#   mean_gap = (WTP[1] − WTP[n]) / (n−1)      (own average WTP step; precomputed at init)
#
#   gap > mean_gap  →  λ *= (1−β)             [inelastic: solidly inframarginal → exploit]
#   gap < mean_gap  →  λ += (1−λ)*α           [elastic: near knife-edge → explore]
#   k=0 or k=n_total  →  hold                 [no next unit to compare]
#
#   Why this equilibrates at λ_eq > 0:
#   At WGE with k_i = n_active*, gap = WTP[n_active*] − WTP[n_active*+1] (the p* step).
#   For uniform-random WTP, this gap is near the median gap, so heat/cool fire roughly
#   equally often → λ stabilises at an intermediate value, not 0.
#   Unlike the rationing signal (only fires when D>S), the elasticity signal fires every
#   turn k_i > 0, giving a balanced equilibrating force independent of market state.
#
#   Crucially: the signal comes only from own WTP array and own k_i.
#   It cannot be corrupted by other agents' draw types (contrast §18 ZIT heat signal).
#
# Init: p̂_i = mean(own WTP), λ_i = 1.0 (pure BiGeo), prev_k_i = 0

using Random
using Statistics
using Printf
using DiscreteMarket

const ETA   = 0.20
const ALPHA = 0.05
const BETA  = 0.35
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

struct ElastConsumer
    wtp      :: Vector{Float64}   # sorted descending
    mean_gap :: Float64           # (wtp[1] - wtp[end]) / (n-1); precomputed
end

struct ElastFirm
    wtac     :: Vector{Float64}   # sorted ascending
    mean_gap :: Float64           # (wtac[end] - wtac[1]) / (n-1); precomputed
end

function ElastConsumer(d::ConsumerDemand)
    w = sort(Float64.(d.wtp), rev=true)
    mg = length(w) > 1 ? (w[1] - w[end]) / (length(w) - 1) : 0.0
    ElastConsumer(w, mg)
end

function ElastFirm(f::FirmSupply)
    c = sort(Float64.(f.wtac))
    mg = length(c) > 1 ? (c[end] - c[1]) / (length(c) - 1) : 0.0
    ElastFirm(c, mg)
end

@inline n_active_c(c::ElastConsumer, p::Float64) = count(v -> v >= p, c.wtp)
@inline n_active_f(f::ElastFirm,     p::Float64) = count(v -> v <= p, f.wtac)

# ── BiGeo draws ───────────────────────────────────────────────────────────────

@inline function draw_bigeo_c(c::ElastConsumer, p_hat::Float64, na::Int,
                               rng::AbstractRNG) :: Int
    na == 0 && return 0
    log_lo = log(max(p_hat, 1e-10))
    log_hi = log(c.wtp[1])
    log_hi <= log_lo + 1e-12 && return na
    pv = exp(log_lo + rand(rng) * (log_hi - log_lo))
    searchsortedlast(c.wtp, pv, rev=true)
end

@inline function draw_bigeo_f(f::ElastFirm, p_hat::Float64, na::Int,
                               rng::AbstractRNG) :: Int
    na == 0 && return 0
    log_lo = log(f.wtac[1])
    log_hi = log(max(p_hat, f.wtac[1] + 1e-10))
    log_hi <= log_lo + 1e-12 && return na
    pv = exp(log_lo + rand(rng) * (log_hi - log_lo))
    searchsortedlast(f.wtac, pv)
end

# ── Mixed draw: BiGeo with prob λ, ratchet with prob 1-λ ─────────────────────
# Returns (d, from_bigeo): demand and draw type.

@inline function draw_c(c::ElastConsumer, p_hat::Float64, λ::Float64,
                         na::Int, prev_k::Int, rng::AbstractRNG) :: Tuple{Int,Bool}
    na == 0 && return (0, false)
    rand(rng) < λ ? (draw_bigeo_c(c, p_hat, na, rng), true) :
                    (min(na, prev_k + 1),               false)
end

@inline function draw_f(f::ElastFirm, p_hat::Float64, λ::Float64,
                         na::Int, prev_k::Int, rng::AbstractRNG) :: Tuple{Int,Bool}
    na == 0 && return (0, false)
    rand(rng) < λ ? (draw_bigeo_f(f, p_hat, na, rng), true) :
                    (min(na, prev_k + 1),               false)
end

# ── p̂ update (BiGeo draws only) ──────────────────────────────────────────────
# Only BiGeo turns update p̂ — ratchet draws (both rationing and at-cap) are excluded.
# This preserves §17's p̂ dynamics exactly: BiGeo provides balanced raise/lower signals
# while the ratchet purely exploits the current p̂ without feeding back noise.

@inline function update_p_c(p::Float64, wtp::Vector{Float64},
                              d::Int, k::Int, na::Int, from_bigeo::Bool) :: Float64
    !from_bigeo && return p
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
    !from_bigeo && return p
    if k < d
        pt = k > 0 ? wtac[k] : wtac[1]
        return (1.0 - ETA) * p + ETA * pt
    elseif d == na && na < length(wtac)
        return (1.0 - ETA) * p + ETA * wtac[na + 1]
    end
    return p
end

# ── λ update (local elasticity gap) ──────────────────────────────────────────

@inline function update_λ_c(λ::Float64, wtp::Vector{Float64},
                              k::Int, mean_gap::Float64) :: Float64
    (k == 0 || k == length(wtp)) && return λ
    gap = wtp[k] - wtp[k + 1]
    gap > mean_gap && return λ * (1.0 - BETA)
    gap < mean_gap && return λ + (1.0 - λ) * ALPHA
    return λ
end

@inline function update_λ_f(λ::Float64, wtac::Vector{Float64},
                              k::Int, mean_gap::Float64) :: Float64
    (k == 0 || k == length(wtac)) && return λ
    gap = wtac[k + 1] - wtac[k]
    gap > mean_gap && return λ * (1.0 - BETA)
    gap < mean_gap && return λ + (1.0 - λ) * ALPHA
    return λ
end

# ── Simulate ──────────────────────────────────────────────────────────────────

function simulate_local_elast(consumers::Vector{ElastConsumer},
                               firms::Vector{ElastFirm},
                               m::GoodMarket, s_wge::Float64, q_star::Int;
                               seed::Int=0)

    rng  = MersenneTwister(seed)
    nc, nf = length(consumers), length(firms)

    p_hat_c  = [mean(c.wtp)  for c in consumers]
    p_hat_f  = [mean(f.wtac) for f in firms]
    lambda_c = ones(Float64, nc)
    lambda_f = ones(Float64, nf)
    prev_k_c = zeros(Int, nc)
    prev_k_f = zeros(Int, nf)

    effs   = Vector{Float64}(undef, T_TURNS)
    q_frac = Vector{Float64}(undef, T_TURNS)
    λ_c_t  = Vector{Float64}(undef, T_TURNS)
    λ_f_t  = Vector{Float64}(undef, T_TURNS)

    for t in 1:T_TURNS
        na_c = [n_active_c(consumers[i], p_hat_c[i]) for i in 1:nc]
        na_f = [n_active_f(firms[j],     p_hat_f[j]) for j in 1:nf]

        drc = [draw_c(consumers[i], p_hat_c[i], lambda_c[i], na_c[i], prev_k_c[i], rng) for i in 1:nc]
        drf = [draw_f(firms[j],     p_hat_f[j], lambda_f[j], na_f[j], prev_k_f[j], rng) for j in 1:nf]

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
            d, bg = drc[i]; k = ks_c[i]; na = na_c[i]
            p_hat_c[i]  = update_p_c(p_hat_c[i], consumers[i].wtp, d, k, na, bg)
            lambda_c[i] = update_λ_c(lambda_c[i], consumers[i].wtp, k, consumers[i].mean_gap)
            !bg && (prev_k_c[i] = k)   # ratchet owns prev_k; BiGeo leaves it unchanged
        end
        for j in 1:nf
            d, bg = drf[j]; k = ks_f[j]; na = na_f[j]
            p_hat_f[j]  = update_p_f(p_hat_f[j], firms[j].wtac, d, k, na, bg)
            lambda_f[j] = update_λ_f(lambda_f[j], firms[j].wtac, k, firms[j].mean_gap)
            !bg && (prev_k_f[j] = k)
        end
    end

    effs, q_frac, λ_c_t, λ_f_t
end

# ── Baseline ──────────────────────────────────────────────────────────────────

function baseline_std(consumers, firms, m, s_wge, rng)
    nc, nf = length(consumers), length(firms)
    mean(begin
        dc = sum(rand(rng, 0:length(consumers[i].wtp)) for i in 1:nc; init=0)
        ds = sum(rand(rng, 0:length(firms[j].wtac))    for j in 1:nf; init=0)
        surplus_f64(m, min(dc, ds)) / s_wge
    end for _ in 1:T_TURNS)
end

# ── Run across markets ────────────────────────────────────────────────────────

function run_local_elast(; nc, nf, n_markets, seed_base, label,
                          max_units=5, report_turns=[1,2,3,5,10,20,50,100,200,300])

    sum_eff   = zeros(T_TURNS); sum_eff2  = zeros(T_TURNS)
    sum_qfrac = zeros(T_TURNS)
    sum_λc    = zeros(T_TURNS); sum_λf    = zeros(T_TURNS)
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

        consumers = [ElastConsumer(d) for d in m.consumers]
        firms     = [ElastFirm(f)     for f in m.firms]

        effs, qf, λc, λf = simulate_local_elast(consumers, firms, m, s_wge, q_star;
                                                  seed=seed + seed_base + 500_000)
        for t in 1:T_TURNS
            sum_eff[t]   += effs[t]; sum_eff2[t]  += effs[t]^2
            sum_qfrac[t] += qf[t]
            sum_λc[t]    += λc[t];   sum_λf[t]    += λf[t]
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
        eff = sum_eff[t] / n; qf  = sum_qfrac[t] / n
        sd  = sqrt(max(0.0, sum_eff2[t]/n - eff^2))
        lc  = sum_λc[t] / n; lf  = sum_λf[t] / n
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

println("§19 LOCAL ELASTICITY ANNEALING — BiGeo + ratchet, elasticity-gap λ signal")
println("  Agents observe ONLY: own WTP/WTAC + own k_i")
println("  State: (p̂_i, λ_i, prev_k_i) per agent")
println()
println("  Quantity:  prob λ   → BiGeo over {WTP ≥ p̂_i}           [explore]")
println("             prob 1-λ → ratchet = min(na(p̂_i), prev_k+1)  [exploit]")
println("  Price:     rationed  → p̂ += η*(WTP[k] - p̂)              [raise]")
println("             at cap    → p̂ += η*(WTP[na+1] - p̂)           [lower]")
println("  Lambda:    gap > mean_gap → λ *= (1-β)                   [cool: inelastic]")
println("             gap < mean_gap → λ += (1-λ)*α                 [heat: elastic]")
println("             gap = WTP[k] - WTP[k+1];  mean_gap precomputed from own WTP")
println("  Init:      p̂ = mean(WTP),  λ = 1.0,  prev_k = 0")
println()

RT = [1, 2, 3, 5, 10, 20, 50, 100, 200, 300]

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

        consumers = [ElastConsumer(d) for d in m.consumers]
        firms     = [ElastFirm(f)     for f in m.firms]

        effs, qf, λc, λf = simulate_local_elast(consumers, firms, m, s_wge, q_star;
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

run_local_elast(nc=50, nf=50, n_markets=300, seed_base=90_000,
                label="100 agents (50C + 50F)  n=300 markets",
                report_turns=RT)

run_local_elast(nc=250, nf=250, n_markets=150, seed_base=91_000,
                label="500 agents (250C + 250F)  n=150 markets",
                report_turns=RT)
