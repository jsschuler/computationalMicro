# Decentralized ZIT — adaptive price beliefs updated from OWN transaction price.
#
# Each agent observes ONLY:
#   - Own WTP / WTAC array (private, known before trading)
#   - Own realized trade k_i (local, per period)
#   - Own transaction price p_t (local: derived from own WTP/WTAC at own realized k_i)
#
# p_t is NOT a global signal. In any real market agents observe what they paid/received.
# They do NOT observe aggregate D, S, q*, or other agents' actions.
#
# Price belief update — THREE cases per period:
#
#   1. Rationed (k < d):
#        Consumer: p_t = WTP[k]         (marginal value received → price is ≥ this)
#                  p̂_i ← EMA(p̂, p_t)   [raise toward WTP[k]]
#        Firm:     p_t = WTAC[k]        (marginal cost supplied → price is ≤ this)
#                  p̂_j ← EMA(p̂, p_t)   [lower toward WTAC[k]]
#
#   2. At ratchet cap, not rationed (k = d = n_active, n_active < n_total):
#        Consumer: p_t = WTP[n_active+1] (first excluded unit → price may be ≤ this)
#                  p̂_i ← EMA(p̂, p_t)   [lower toward WTP[n_active+1]]
#        Firm:     p_t = WTAC[n_active+1] (first excluded unit → price may be ≥ this)
#                  p̂_j ← EMA(p̂, p_t)   [raise toward WTAC[n_active+1]]
#
#   3. Ratcheting (k = d < n_active) OR at full capacity (n_active = n_total):
#        No update — no informative price signal available.
#
# Convergence: p̂_i oscillates between WTP[n_active*] and WTP[n_active*+1], centred near p*.
#              Amplitude = WTP gap straddling p* → 0 as WTP distribution becomes dense.
#
# Quantity: ratchet  d_i(t) = min(count(WTP ≥ p̂_i), k_i(t-1) + 1)
# Init:     p̂_i = mean(own WTP),   p̂_j = mean(own WTAC)   [neutral prior]
#
# Note on Agda formalization (see §16):
#   The false-equilibrium trap of §16 required only {rationing, own WTP/WTAC} as signals.
#   This mechanism adds "own transaction price" = the value of the own marginal traded unit.
#   The minimal information requirement for WGE convergence is therefore:
#     LocalInfo = (own_WTP, own_WTAC, own_k_i, own_p_t)
#   where p_t is DERIVED from own_WTP[k_i] — no additional global signal needed.

using Random
using Statistics
using Printf
using DiscreteMarket

const ETA_RAT = 0.20   # EMA weight when rationed
const ETA_CAP = 0.20   # EMA weight when at cap, not rationed (symmetric)
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

struct DecentConsumer
    wtp :: Vector{Float64}   # all valuations, sorted descending
end

struct DecentFirm
    wtac :: Vector{Float64}  # all costs, sorted ascending
end

DecentConsumer(d::ConsumerDemand) = DecentConsumer(sort(Float64.(d.wtp), rev=true))
DecentFirm(f::FirmSupply)         = DecentFirm(sort(Float64.(f.wtac)))

@inline n_active_c(c::DecentConsumer, p::Float64) = count(v -> v >= p, c.wtp)
@inline n_active_f(f::DecentFirm,     p::Float64) = count(v -> v <= p, f.wtac)

# ── Price belief update ───────────────────────────────────────────────────────

# Consumer:
#   Rationed (k<d):             pull p̂ UP   toward WTP[k]       (price is ≥ WTP[k])
#   At cap, not rationed:       pull p̂ DOWN toward WTP[na+1]    (next excluded unit)
#   Ratcheting / at n_total cap: no update
@inline function update_p_c(p::Float64, wtp::Vector{Float64},
                              d::Int, k::Int, na::Int) :: Float64
    if k < d
        # rationed: price ≥ WTP[k] → raise belief (strong correction)
        pt = k > 0 ? wtp[k] : wtp[1]
        return (1.0 - ETA_RAT) * p + ETA_RAT * pt
    elseif d == na && na < length(wtp)
        # at cap, not rationed: price ≤ WTP[na+1] → lower belief (gentle exploration)
        pt = wtp[na + 1]
        return (1.0 - ETA_CAP) * p + ETA_CAP * pt
    end
    return p   # ratcheting or at full capacity → no signal
end

# Firm (symmetric, reversed directions):
@inline function update_p_f(p::Float64, wtac::Vector{Float64},
                              d::Int, k::Int, na::Int) :: Float64
    if k < d
        pt = k > 0 ? wtac[k] : wtac[1]
        return (1.0 - ETA_RAT) * p + ETA_RAT * pt
    elseif d == na && na < length(wtac)
        pt = wtac[na + 1]
        return (1.0 - ETA_CAP) * p + ETA_CAP * pt
    end
    return p
end

# ── Simulate ──────────────────────────────────────────────────────────────────

function simulate_decent_price(consumers::Vector{DecentConsumer},
                                firms::Vector{DecentFirm},
                                m::GoodMarket, s_wge::Float64, q_star::Int) :: Tuple{Vector{Float64},Vector{Float64}}

    nc, nf = length(consumers), length(firms)

    # Neutral prior: mean of own valuations
    p_hat_c  = [mean(c.wtp) for c in consumers]
    p_hat_f  = [mean(f.wtac) for f in firms]
    prev_k_c = zeros(Int, nc)
    prev_k_f = zeros(Int, nf)

    effs   = Vector{Float64}(undef, T_TURNS)
    q_frac = Vector{Float64}(undef, T_TURNS)

    for t in 1:T_TURNS
        na_c = [n_active_c(consumers[i], p_hat_c[i]) for i in 1:nc]
        na_f = [n_active_f(firms[j],     p_hat_f[j]) for j in 1:nf]

        ds_c = [min(na_c[i], prev_k_c[i] + 1) for i in 1:nc]
        ds_f = [min(na_f[j], prev_k_f[j] + 1) for j in 1:nf]

        D = sum(ds_c)
        S = sum(ds_f)
        q = min(D, S)

        effs[t]   = surplus_f64(m, q) / s_wge
        q_frac[t] = q_star > 0 ? q / q_star : 1.0

        ks_c = [realized_c(ds_c[i], D, S) for i in 1:nc]
        ks_f = [realized_f(ds_f[j], D, S) for j in 1:nf]

        for i in 1:nc
            p_hat_c[i]  = update_p_c(p_hat_c[i], consumers[i].wtp, ds_c[i], ks_c[i], na_c[i])
            prev_k_c[i] = ks_c[i]
        end
        for j in 1:nf
            p_hat_f[j]  = update_p_f(p_hat_f[j], firms[j].wtac, ds_f[j], ks_f[j], na_f[j])
            prev_k_f[j] = ks_f[j]
        end
    end

    effs, q_frac
end

# ── Baseline: standard ZIT ────────────────────────────────────────────────────

@inline function baseline_std(consumers, firms, m, s_wge, T, rng)
    nc, nf = length(consumers), length(firms)
    mean(begin
        dc = sum(rand(rng, 0:length(consumers[i].wtp)) for i in 1:nc; init=0)
        ds = sum(rand(rng, 0:length(firms[j].wtac))    for j in 1:nf; init=0)
        surplus_f64(m, min(dc, ds)) / s_wge
    end for _ in 1:T)
end

# ── Run across markets ────────────────────────────────────────────────────────

function run_decent_price(; nc, nf, n_markets, seed_base, label,
                           max_units=5, report_turns=[1,2,3,5,10,20,50,100,200,300])

    sum_eff   = zeros(T_TURNS)
    sum_eff2  = zeros(T_TURNS)
    sum_qfrac = zeros(T_TURNS)
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

        consumers = [DecentConsumer(d) for d in m.consumers]
        firms     = [DecentFirm(f)     for f in m.firms]

        effs, qf = simulate_decent_price(consumers, firms, m, s_wge, q_star)
        for t in 1:T_TURNS
            sum_eff[t]   += effs[t]
            sum_eff2[t]  += effs[t]^2
            sum_qfrac[t] += qf[t]
        end
        n_analyzed += 1

        rng2 = MersenneTwister(seed + seed_base + 99_000)
        push!(std_effs, baseline_std(consumers, firms, m, s_wge, T_TURNS, rng2))
    end

    n = n_analyzed
    println("=" ^ 70)
    println(label)
    @printf("  n=%d  T=%d  ETA_RAT=%.2f, ETA_CAP=%.2f\n\n", n, T_TURNS, ETA_RAT, ETA_CAP)
    @printf("  Baseline — Static Std:   mean eff = %.4f\n", mean(std_effs))
    @printf("  Target   — WGE:          mean eff = 1.0000\n")
    println()
    println("  Turn │  q/q*    Eff_mean  Eff_std  vs Std   vs WGE")
    println("  ─────┼──────────────────────────────────────────────")
    for t in report_turns
        t > T_TURNS && continue
        eff = sum_eff[t] / n
        qf  = sum_qfrac[t] / n
        sd  = sqrt(max(0.0, sum_eff2[t]/n - eff^2))
        @printf("  t=%3d │  %.4f  %.4f    %.4f   %+.4f  %+.4f\n",
                t, qf, eff, sd, eff - mean(std_effs), eff - 1.0)
    end
    println()
    final_eff = sum_eff[T_TURNS]  / n
    final_sd  = sqrt(max(0.0, sum_eff2[T_TURNS]/n - final_eff^2))
    final_qf  = sum_qfrac[T_TURNS] / n
    @printf("  Final (t=%d):\n", T_TURNS)
    @printf("    q/q* = %.4f\n", final_qf)
    @printf("    Efficiency: mean=%.4f  std=%.4f  vs WGE=%+.4f\n",
            final_eff, final_sd, final_eff - 1.0)
    println()

    (sum_eff ./ n, sum_qfrac ./ n, n)
end

# ── Main ──────────────────────────────────────────────────────────────────────

println("DECENTRALIZED ZIT — transaction-price adaptive beliefs (corrected)")
println("  Local signals: own WTP/WTAC + own k_i + own p_t derived from WTP[k_i]")
println()
println("  Update rules:")
println("    Rationed (k<d):          p̂_c += η*(WTP[k] - p̂_c)       [raise: price ≥ WTP[k]]")
println("    At cap, not ration:      p̂_c += η*(WTP[na+1] - p̂_c)    [lower: price ≤ WTP[na+1]]")
println("    Ratcheting / full cap:   no update                        [no informative signal]")
println("    Firm rationed:           p̂_f += η*(WTAC[k] - p̂_f)      [lower: price ≤ WTAC[k]]")
println("    Firm at cap, not ration: p̂_f += η*(WTAC[na+1] - p̂_f)  [raise: price ≥ WTAC[na+1]]")
println("  Convergence: p̂_i oscillates between WTP[na*] and WTP[na*+1], centred near p*")
println("  Init: p̂_i = mean(own WTP), p̂_j = mean(own WTAC)")
println()

RT = [1, 2, 3, 5, 10, 20, 50, 100, 200, 300]

# Small / mixed markets
let
    SEED_BASE  = 60_000
    sum_eff    = zeros(T_TURNS); sum_eff2  = zeros(T_TURNS)
    sum_qfrac  = zeros(T_TURNS)
    n          = 0
    std_effs   = Float64[]

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

        consumers = [DecentConsumer(d) for d in m.consumers]
        firms     = [DecentFirm(f)     for f in m.firms]

        effs, qf = simulate_decent_price(consumers, firms, m, s_wge, q_star)
        for t in 1:T_TURNS
            sum_eff[t]  += effs[t]; sum_eff2[t]  += effs[t]^2
            sum_qfrac[t] += qf[t]
        end
        n += 1

        rng2 = MersenneTwister(seed + SEED_BASE + 99_000)
        push!(std_effs, baseline_std(consumers, firms, m, s_wge, T_TURNS, rng2))
    end

    println("=" ^ 70)
    println("Small/mixed markets  (1–8C, 1–8F, 1–6 max_units)  n=$n / 1000")
    @printf("  T=%d  ETA_RAT=%.2f, ETA_CAP=%.2f\n\n", T_TURNS, ETA_RAT, ETA_CAP)
    @printf("  Baseline — Static Std:   mean eff = %.4f\n", mean(std_effs))
    println()
    println("  Turn │  q/q*    Eff_mean  Eff_std  vs Std   vs WGE")
    println("  ─────┼──────────────────────────────────────────────")
    for t in RT
        t > T_TURNS && continue
        eff = sum_eff[t] / n; qf = sum_qfrac[t] / n
        sd  = sqrt(max(0.0, sum_eff2[t]/n - eff^2))
        @printf("  t=%3d │  %.4f  %.4f    %.4f   %+.4f  %+.4f\n",
                t, qf, eff, sd, eff - mean(std_effs), eff - 1.0)
    end
    println()
    fe = sum_eff[T_TURNS]/n; fq = sum_qfrac[T_TURNS]/n
    @printf("  Final (t=%d): q/q*=%.4f  eff=%.4f  std=%.4f  vs WGE=%+.4f\n\n",
            T_TURNS, fq, fe, sqrt(max(0.0,sum_eff2[T_TURNS]/n - fe^2)), fe - 1.0)
end

run_decent_price(nc=50, nf=50, n_markets=300, seed_base=70_000,
                 label="100 agents (50C + 50F)  n=300 markets",
                 report_turns=RT)

run_decent_price(nc=250, nf=250, n_markets=150, seed_base=71_000,
                 label="500 agents (250C + 250F)  n=150 markets",
                 report_turns=RT)
