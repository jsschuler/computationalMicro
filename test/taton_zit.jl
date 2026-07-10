# Rationing-guided quantity tatonnement — convergence to WGE without global signals.
#
# Each agent maintains prev_realized (last period's realized trade, starts at 0).
# Each period, demand/supply = min(n_active_i, prev_realized + 1).
#
# This is a RATCHET: always probe one more unit than the market last allowed,
# capped at n_active_i so agents never demand unprofitable units.
#
# The ONLY local information used:
#   - Own realized trade k_i each period (private feedback from pro-rata allocation)
#   - Own WTP / WTAC array and p* (private knowledge, used to compute n_active_i)
#   - No global D, S, or q* ever observed
#
# Why it converges:
#   - prev_realized ≤ n_active_i always (cap prevents overdemand)
#   - d_i = prev_realized + 1 ≤ n_active_i → agents always want one more unit if available
#   - When D ≤ S: k_i = d_i → prev_realized = d_i → next d = d_i + 1 (increment)
#   - When D > S: k_i < d_i → prev_realized = k_i < d_i → next d = k_i + 1 ≤ d_i (hold or shrink)
#   - Fixed point: prev_realized = n_active_i, d_i = n_active_i + 1... wait, no:
#     at WGE D=S: k_i = d_i = n_active_i. prev_realized = n_active_i.
#     next d = min(n_active_i, n_active_i + 1) = n_active_i. HOLD. ✓
#
# Stochastic comparison: also runs pure BiGeo and pure standard ZIT baselines.
#
# Run with:
#   julia --project=. test/taton_zit.jl

using Random
using Statistics
using Printf
using DiscreteMarket

const T_TURNS   = 200
const MAX_UNITS = 5

# ── Market-level surplus ──────────────────────────────────────────────────────

function surplus_f64(m::GoodMarket, q::Int) :: Float64
    q <= 0 && return 0.0
    all_wtp  = sort([Float64(v) for d in m.consumers for v in d.wtp], rev=true)
    all_wtac = sort([Float64(c) for f in m.firms     for c in f.wtac])
    n = min(q, length(all_wtp), length(all_wtac))
    n == 0 && return 0.0
    sum(all_wtp[1:n]) - sum(all_wtac[1:n])
end

# ── Precomputed per-agent data ─────────────────────────────────────────────────

struct TatonConsumer
    n_active :: Int              # WGE-optimal demand = count(wtp ≥ p*)
    n_total  :: Int              # upper bound on draw
    log_p    :: Float64
    log_hi   :: Float64
    active   :: Vector{Float64}  # wtp ≥ p*, sorted descending
end

struct TatonFirm
    n_active :: Int
    n_total  :: Int
    log_lo   :: Float64
    log_p    :: Float64
    active   :: Vector{Float64}  # wtac ≤ p*, sorted ascending
end

function TatonConsumer(d::ConsumerDemand, p_star::Price)
    lp  = log(Float64(p_star))
    act = sort([Float64(v) for v in d.wtp if v >= p_star], rev=true)
    lhi = isempty(act) ? lp : log(act[1])
    TatonConsumer(length(act), length(d.wtp), lp, lhi, act)
end

function TatonFirm(f::FirmSupply, p_star::Price)
    lp  = log(Float64(p_star))
    act = sort([Float64(c) for c in f.wtac if c <= p_star])
    llo = isempty(act) ? lp : log(act[1])
    TatonFirm(length(act), length(f.wtac), llo, lp, act)
end

# ── Pro-rata realized trade ───────────────────────────────────────────────────

@inline realized_c(draw, D, S) = D == 0 ? 0 : (D <= S ? draw : floor(Int, draw * S / D))
@inline realized_f(draw, D, S) = S == 0 ? 0 : (S <= D ? draw : floor(Int, draw * D / S))

# ── Ratchet demand/supply from realized ──────────────────────────────────────
# Next draw = min(n_active, prev_realized + 1)
# Agent probes one more unit than the market last allowed, capped at WGE-optimal.

@inline ratchet(prev_k::Int, n_active::Int) :: Int = min(n_active, prev_k + 1)

# ── BiGeo draw (for baseline comparison) ─────────────────────────────────────

@inline function draw_bgip_consumer(c::TatonConsumer, rng::AbstractRNG) :: Int
    isempty(c.active) && return 0
    c.log_hi ≈ c.log_p && return length(c.active)
    pv = exp(c.log_p + rand(rng) * (c.log_hi - c.log_p))
    searchsortedlast(c.active, pv, rev=true)
end

@inline function draw_bgip_firm(f::TatonFirm, rng::AbstractRNG) :: Int
    isempty(f.active) && return 0
    f.log_p ≈ f.log_lo && return length(f.active)
    pv = exp(f.log_lo + rand(rng) * (f.log_p - f.log_lo))
    searchsortedlast(f.active, pv)
end

# ── Simulate tatonnement ──────────────────────────────────────────────────────

function simulate_taton(consumers::Vector{TatonConsumer},
                         firms::Vector{TatonFirm},
                         m::GoodMarket, s_wge::Float64, q_star::Int;
                         seed::Int=0) :: Tuple{Vector{Float64},Vector{Float64}}

    nc, nf = length(consumers), length(firms)

    # Ratchet state: each agent's last realized trade (initialized to 0)
    # Next demand = min(n_active_i, prev_k + 1)
    prev_k_c = zeros(Int, nc)
    prev_k_f = zeros(Int, nf)

    effs   = Vector{Float64}(undef, T_TURNS)
    q_frac = Vector{Float64}(undef, T_TURNS)

    for t in 1:T_TURNS
        # Demands and supplies this period: one more than last period's realized
        D = sum(ratchet(prev_k_c[i], consumers[i].n_active) for i in 1:nc)
        S = sum(ratchet(prev_k_f[j], firms[j].n_active)     for j in 1:nf)
        q = min(D, S)

        effs[t]   = surplus_f64(m, q) / s_wge
        q_frac[t] = q_star > 0 ? q / q_star : 1.0

        # Observe realized and store as prev_k for next period
        for i in 1:nc
            d = ratchet(prev_k_c[i], consumers[i].n_active)
            prev_k_c[i] = realized_c(d, D, S)
        end
        for j in 1:nf
            d = ratchet(prev_k_f[j], firms[j].n_active)
            prev_k_f[j] = realized_f(d, D, S)
        end
    end

    effs, q_frac
end

# ── Run across markets ────────────────────────────────────────────────────────

function run_taton(; nc, nf, n_markets, seed_base, label, max_units=MAX_UNITS,
                    report_turns=[1,2,3,5,10,20,30,50,100,150,200])

    sum_eff   = zeros(T_TURNS)
    sum_eff2  = zeros(T_TURNS)
    sum_qfrac = zeros(T_TURNS)
    n_analyzed = 0
    static_std   = Float64[]
    static_bgip  = Float64[]

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

        consumers = [TatonConsumer(d, p_star) for d in m.consumers]
        firms     = [TatonFirm(f, p_star)     for f in m.firms]

        effs, qf = simulate_taton(consumers, firms, m, s_wge, q_star; seed=seed)
        for t in 1:T_TURNS
            sum_eff[t]   += effs[t]
            sum_eff2[t]  += effs[t]^2
            sum_qfrac[t] += qf[t]
        end
        n_analyzed += 1

        # Baselines
        rng2 = MersenneTwister(seed + seed_base + 99_000)
        push!(static_std, mean(begin
            dc = sum(rand(rng2, 0:consumers[i].n_total) for i in 1:nc; init=0)
            ds = sum(rand(rng2, 0:firms[j].n_total)     for j in 1:nf; init=0)
            surplus_f64(m, min(dc, ds)) / s_wge
        end for _ in 1:T_TURNS))

        rng3 = MersenneTwister(seed + seed_base + 199_000)
        push!(static_bgip, mean(begin
            dc = sum(draw_bgip_consumer(consumers[i], rng3) for i in 1:nc; init=0)
            ds = sum(draw_bgip_firm(firms[j], rng3)         for j in 1:nf; init=0)
            surplus_f64(m, min(dc, ds)) / s_wge
        end for _ in 1:T_TURNS))
    end

    n = n_analyzed
    println("=" ^ 70)
    println(label)
    @printf("  n=%d  T=%d  start=d_i=0 for all agents\n\n", n, T_TURNS)
    @printf("  Baseline — Static Std:   mean eff = %.4f\n", mean(static_std))
    @printf("  Baseline — Static BiGeo: mean eff = %.4f\n", mean(static_bgip))
    @printf("  Target   — WGE:          mean eff = 1.0000\n")
    println()

    println("  Turn │  q/q*    Eff_mean  Eff_std  vs Std   vs BiGeo  vs WGE")
    println("  ─────┼────────────────────────────────────────────────────────────")
    for t in report_turns
        t > T_TURNS && continue
        eff  = sum_eff[t]   / n
        qf   = sum_qfrac[t] / n
        var  = sum_eff2[t]/n - eff^2
        sd   = var > 0 ? sqrt(var) : 0.0
        @printf("  t=%3d │  %.4f  %.4f    %.4f   %+.4f  %+.4f   %+.4f\n",
                t, qf, eff, sd,
                eff - mean(static_std), eff - mean(static_bgip), eff - 1.0)
    end
    println()

    final_eff = sum_eff[T_TURNS]  / n
    final_std = sqrt(max(0.0, sum_eff2[T_TURNS]/n - final_eff^2))
    final_qf  = sum_qfrac[T_TURNS] / n
    @printf("  Final (t=%d):\n", T_TURNS)
    @printf("    q/q* = %.4f  (%.1f%% of WGE quantity)\n", final_qf, 100*final_qf)
    @printf("    Efficiency: mean=%.4f  std=%.4f\n", final_eff, final_std)
    @printf("    vs WGE:  %+.4f\n", final_eff - 1.0)
    @printf("    vs Std:  %+.4f\n", final_eff - mean(static_std))
    @printf("    vs BiGeo:%+.4f\n", final_eff - mean(static_bgip))
    println()

    (sum_eff ./ n, sum_qfrac ./ n, n)
end

# ── Main ──────────────────────────────────────────────────────────────────────

println("RATIONING-GUIDED TATONNEMENT — Convergence to WGE without global signals")
println("  Local signals only: own realized trade k_i, own WTP/WTAC, own p*")
println("  No D, S, q*, or aggregate efficiency observed")
println()
println("  Ratchet rule:  d_i(t) = min(n_active_i,  k_i(t-1) + 1)")
println("    Never demands unprofitable units (cap at n_active_i)")
println("    Always probes one more than the market last allowed")
println("    Fixed point: k_i = n_active_i → d_i = n_active_i → WGE")
println("    Initialized: prev_k = 0 for all agents")
println()

RT = [1, 2, 3, 5, 10, 20, 30, 50, 100, 150, 200]

# Small / mixed markets
let
    SEED_BASE = 20_000
    sum_eff   = zeros(T_TURNS); sum_eff2  = zeros(T_TURNS)
    sum_qfrac = zeros(T_TURNS)
    n = 0
    static_std  = Float64[]; static_bgip = Float64[]

    for seed in 1:1_000
        rng   = MersenneTwister(seed + SEED_BASE)
        n_c   = rand(rng, 1:8); n_f = rand(rng, 1:8); max_u = rand(rng, 1:6)
        m, r  = generate_good_market(rng;
            good=1, n_consumers=n_c, n_firms=n_f, max_units=max_u, Q=100)
        r.cleared || continue
        p_star = r.price
        q_star = aggregate_demand(m, p_star)
        s_wge  = surplus_f64(m, q_star)
        s_wge > 0 || continue

        consumers = [TatonConsumer(d, p_star) for d in m.consumers]
        firms     = [TatonFirm(f, p_star)     for f in m.firms]

        effs, qf = simulate_taton(consumers, firms, m, s_wge, q_star; seed=seed)
        for t in 1:T_TURNS
            sum_eff[t]   += effs[t]; sum_eff2[t]  += effs[t]^2
            sum_qfrac[t] += qf[t]
        end
        n += 1

        rng2 = MersenneTwister(seed + SEED_BASE + 99_000)
        push!(static_std, mean(begin
            dc = sum(rand(rng2, 0:consumers[i].n_total) for i in 1:n_c; init=0)
            ds = sum(rand(rng2, 0:firms[j].n_total)     for j in 1:n_f; init=0)
            surplus_f64(m, min(dc, ds)) / s_wge end for _ in 1:T_TURNS))
        rng3 = MersenneTwister(seed + SEED_BASE + 199_000)
        push!(static_bgip, mean(begin
            dc = sum(draw_bgip_consumer(consumers[i], rng3) for i in 1:n_c; init=0)
            ds = sum(draw_bgip_firm(firms[j], rng3)         for j in 1:n_f; init=0)
            surplus_f64(m, min(dc, ds)) / s_wge end for _ in 1:T_TURNS))
    end

    println("=" ^ 70)
    println("Small/mixed markets  (1–8C, 1–8F, 1–6 max_units)  n=$n / 1000")
    @printf("  T=%d\n\n", T_TURNS)
    @printf("  Baseline — Static Std:   mean eff = %.4f\n", mean(static_std))
    @printf("  Baseline — Static BiGeo: mean eff = %.4f\n", mean(static_bgip))
    println()
    println("  Turn │  q/q*    Eff_mean  Eff_std  vs Std   vs BiGeo  vs WGE")
    println("  ─────┼────────────────────────────────────────────────────────────")
    for t in RT
        t > T_TURNS && continue
        eff = sum_eff[t]  / n;  qf = sum_qfrac[t] / n
        sd  = sqrt(max(0.0, sum_eff2[t]/n - eff^2))
        @printf("  t=%3d │  %.4f  %.4f    %.4f   %+.4f  %+.4f   %+.4f\n",
                t, qf, eff, sd, eff-mean(static_std), eff-mean(static_bgip), eff-1.0)
    end
    println()
    final_eff = sum_eff[T_TURNS]/n;  final_qf = sum_qfrac[T_TURNS]/n
    @printf("  Final (t=%d): q/q*=%.4f  eff=%.4f  std=%.4f  vs WGE=%+.4f\n\n",
            T_TURNS, final_qf, final_eff,
            sqrt(max(0.0, sum_eff2[T_TURNS]/n - final_eff^2)), final_eff-1.0)
end

# 100-agent markets
run_taton(nc=50, nf=50, n_markets=300, seed_base=30_000,
          label="100 agents (50C + 50F)  n=300 markets",
          report_turns=RT)

# 500-agent markets
run_taton(nc=250, nf=250, n_markets=150, seed_base=31_000,
          label="500 agents (250C + 250F)  n=150 markets",
          report_turns=RT)
