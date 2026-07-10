# Adaptive ZIT: agents learn to switch between additive and BiGeo strategies.
#
# Each agent independently tracks the last observed performance under each
# strategy (standard ZIT and BiGeo price-space ZIT).  Each turn:
#   1. Agent plays current strategy, observes realized trade quantity.
#   2. Updates performance memory for that strategy.
#   3. Switches to whichever strategy has the better remembered performance.
#      (NaN = never tried → treated as -Inf, forcing exploration on turn 2.)
#
# Cold start: all agents begin with :std and NaN memory for both strategies.
# Turn 1: play :std, get perf.  BiGeo memory = NaN → -Inf < perf_std → switch.
# Turn 2: all play :bgip.  From turn 3 onward, genuine comparison rules.
#
# Performance signal: pro-rata realized trade quantity.
#   Consumers: if D ≤ S, get full draw; if D > S, get floor(draw * S/D).
#   Firms:     if S ≤ D, get full draw; if S > D, get floor(draw * D/S).
#
# Run with:
#   julia --project=. test/adaptive_zit.jl

using Random
using Statistics
using Printf
using DiscreteMarket

const T_TURNS  = 300
const MAX_UNITS = 5

# ── Surplus ───────────────────────────────────────────────────────────────────

function surplus_f64(m::GoodMarket, q::Int) :: Float64
    q <= 0 && return 0.0
    all_wtp  = sort([Float64(v) for d in m.consumers for v in d.wtp], rev=true)
    all_wtac = sort([Float64(c) for f in m.firms     for c in f.wtac])
    n = min(q, length(all_wtp), length(all_wtac))
    n == 0 && return 0.0
    sum(all_wtp[1:n]) - sum(all_wtac[1:n])
end

# ── Precomputed per-agent data (same structs as bgi_price_large.jl) ───────────

struct AdaptConsumer
    n_total  :: Int              # standard ZIT max_q = len(wtp)
    log_p    :: Float64
    log_hi   :: Float64          # log(max active wtp); = log_p if no spread
    active   :: Vector{Float64}  # wtp values ≥ p*, sorted descending
end

struct AdaptFirm
    n_total  :: Int
    log_lo   :: Float64          # log(min active wtac); = log_p if no spread
    log_p    :: Float64
    active   :: Vector{Float64}  # wtac values ≤ p*, sorted ascending
end

function AdaptConsumer(d::ConsumerDemand, p_star::Price)
    lp  = log(Float64(p_star))
    act = sort([Float64(v) for v in d.wtp if v >= p_star], rev=true)
    lhi = isempty(act) ? lp : log(act[1])
    AdaptConsumer(length(d.wtp), lp, lhi, act)
end

function AdaptFirm(f::FirmSupply, p_star::Price)
    lp  = log(Float64(p_star))
    act = sort([Float64(c) for c in f.wtac if c <= p_star])
    llo = isempty(act) ? lp : log(act[1])
    AdaptFirm(length(f.wtac), llo, lp, act)
end

# ── Single-agent draws ────────────────────────────────────────────────────────

@inline function draw_std_consumer(c::AdaptConsumer, rng::AbstractRNG) :: Int
    rand(rng, 0:c.n_total)
end

@inline function draw_bgip_consumer(c::AdaptConsumer, rng::AbstractRNG) :: Int
    isempty(c.active) && return 0
    c.log_hi ≈ c.log_p && return length(c.active)
    pv = exp(c.log_p + rand(rng) * (c.log_hi - c.log_p))
    searchsortedlast(c.active, pv, rev=true)
end

@inline function draw_std_firm(f::AdaptFirm, rng::AbstractRNG) :: Int
    rand(rng, 0:f.n_total)
end

@inline function draw_bgip_firm(f::AdaptFirm, rng::AbstractRNG) :: Int
    isempty(f.active) && return 0
    f.log_p ≈ f.log_lo && return length(f.active)
    pv = exp(f.log_lo + rand(rng) * (f.log_p - f.log_lo))
    searchsortedlast(f.active, pv)
end

# ── Pro-rata realized trade ───────────────────────────────────────────────────

@inline realized_c(draw, D, S) = D == 0 ? 0 : (D <= S ? draw : floor(Int, draw * S / D))
@inline realized_f(draw, D, S) = S == 0 ? 0 : (S <= D ? draw : floor(Int, draw * D / S))

# ── Next-strategy rule ────────────────────────────────────────────────────────
# NaN memory → -Inf (never tried → always appears worse → forces exploration)

@inline function next_strategy(perf_std::Float64, perf_bgip::Float64) :: Bool
    # Untried strategy gets +Inf (optimistic prior) → always explore it once.
    # Tie → stay std (false).
    ps = isnan(perf_std)  ? +Inf : perf_std
    pb = isnan(perf_bgip) ? +Inf : perf_bgip
    pb > ps
end

# ── Simulate one market for T_TURNS turns ────────────────────────────────────
# Returns (λ_c[t], λ_f[t], eff[t]) vectors over turns.

function simulate_market(consumers::Vector{AdaptConsumer},
                          firms::Vector{AdaptFirm},
                          m::GoodMarket, s_wge::Float64;
                          seed::Int=0) :: Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}

    nc, nf = length(consumers), length(firms)

    # Agent state: strategy (false=std, true=bgip) + performance memory
    use_bgip_c   = falses(nc)
    perf_std_c   = fill(NaN, nc)
    perf_bgip_c  = fill(NaN, nc)

    use_bgip_f   = falses(nf)
    perf_std_f   = fill(NaN, nf)
    perf_bgip_f  = fill(NaN, nf)

    λ_c  = Vector{Float64}(undef, T_TURNS)
    λ_f  = Vector{Float64}(undef, T_TURNS)
    effs = Vector{Float64}(undef, T_TURNS)

    rng = MersenneTwister(seed)

    draws_c = Vector{Int}(undef, nc)
    draws_f = Vector{Int}(undef, nf)

    for t in 1:T_TURNS
        # ── Draw ─────────────────────────────────────────────────────────────
        for i in 1:nc
            draws_c[i] = use_bgip_c[i] ?
                draw_bgip_consumer(consumers[i], rng) :
                draw_std_consumer(consumers[i], rng)
        end
        for j in 1:nf
            draws_f[j] = use_bgip_f[j] ?
                draw_bgip_firm(firms[j], rng) :
                draw_std_firm(firms[j], rng)
        end

        D = sum(draws_c)
        S = sum(draws_f)

        effs[t] = surplus_f64(m, min(D, S)) / s_wge
        λ_c[t]  = mean(use_bgip_c)
        λ_f[t]  = mean(use_bgip_f)

        # ── Update performance memory ─────────────────────────────────────────
        for i in 1:nc
            r = realized_c(draws_c[i], D, S)
            if use_bgip_c[i]
                perf_bgip_c[i] = Float64(r)
            else
                perf_std_c[i]  = Float64(r)
            end
        end
        for j in 1:nf
            r = realized_f(draws_f[j], D, S)
            if use_bgip_f[j]
                perf_bgip_f[j] = Float64(r)
            else
                perf_std_f[j]  = Float64(r)
            end
        end

        # ── Decide next strategy ──────────────────────────────────────────────
        for i in 1:nc
            use_bgip_c[i] = next_strategy(perf_std_c[i], perf_bgip_c[i])
        end
        for j in 1:nf
            use_bgip_f[j] = next_strategy(perf_std_f[j], perf_bgip_f[j])
        end
    end

    λ_c, λ_f, effs
end

# ── Run across many markets, report turn-by-turn averages ────────────────────

function run_adaptive(; nc, nf, n_markets, seed_base, label, max_units=MAX_UNITS,
                       report_turns=[1,2,3,5,10,25,50,100,200,300])

    # Accumulators: sum over markets at each turn
    sum_λ_c  = zeros(T_TURNS)
    sum_λ_f  = zeros(T_TURNS)
    sum_eff  = zeros(T_TURNS)
    sum_eff2 = zeros(T_TURNS)
    n_analyzed = 0

    # Also accumulate static baselines across same markets
    static_eff_std  = Float64[]
    static_eff_bgip = Float64[]

    for seed in 1:n_markets
        rng = MersenneTwister(seed + seed_base)
        m, r = generate_good_market(rng;
            good=1, n_consumers=nc, n_firms=nf,
            max_units=max_units, Q=100)
        r.cleared || continue
        p_star = r.price
        s_wge  = surplus_f64(m, aggregate_demand(m, p_star))
        s_wge > 0 || continue

        consumers = [AdaptConsumer(d, p_star) for d in m.consumers]
        firms     = [AdaptFirm(f, p_star)     for f in m.firms]

        λ_c, λ_f, effs = simulate_market(consumers, firms, m, s_wge; seed=seed)

        for t in 1:T_TURNS
            sum_λ_c[t]  += λ_c[t]
            sum_λ_f[t]  += λ_f[t]
            sum_eff[t]  += effs[t]
            sum_eff2[t] += effs[t]^2
        end
        n_analyzed += 1

        # Static baselines: run one block of T_TURNS turns each, pure strategies
        rng2 = MersenneTwister(seed + seed_base + 99_000)
        s_std = mean(begin
            dc = sum(draw_std_consumer(consumers[i], rng2) for i in 1:nc; init=0)
            ds = sum(draw_std_firm(firms[j], rng2) for j in 1:nf; init=0)
            surplus_f64(m, min(dc, ds)) / s_wge
        end for _ in 1:T_TURNS)

        rng3 = MersenneTwister(seed + seed_base + 199_000)
        s_bgip = mean(begin
            dc = sum(draw_bgip_consumer(consumers[i], rng3) for i in 1:nc; init=0)
            ds = sum(draw_bgip_firm(firms[j], rng3) for j in 1:nf; init=0)
            surplus_f64(m, min(dc, ds)) / s_wge
        end for _ in 1:T_TURNS)

        push!(static_eff_std,  s_std)
        push!(static_eff_bgip, s_bgip)
    end

    n = n_analyzed
    println("=" ^ 66)
    println(label)
    @printf("  n_analyzed=%d / %d  |  T_TURNS=%d\n\n", n, n_markets, T_TURNS)

    # Static baselines
    @printf("  Static λ=0 (pure Std):   mean eff = %.4f\n", mean(static_eff_std))
    @printf("  Static λ=1 (pure BiGeo): mean eff = %.4f\n", mean(static_eff_bgip))
    println()

    println("  Turn │  λ_c    λ_f    Eff_mean  Eff_std  vs Std   vs BiGeo")
    println("  ─────┼────────────────────────────────────────────────────────")

    for t in report_turns
        t > T_TURNS && continue
        λc  = sum_λ_c[t] / n
        λf  = sum_λ_f[t] / n
        eff = sum_eff[t]  / n
        var = sum_eff2[t] / n - eff^2
        sd  = var > 0 ? sqrt(var) : 0.0
        Δstd  = eff - mean(static_eff_std)
        Δbgip = eff - mean(static_eff_bgip)
        @printf("  t=%3d │  %.3f  %.3f  %.4f    %.4f   %+.4f  %+.4f\n",
                t, λc, λf, eff, sd, Δstd, Δbgip)
    end
    println()

    # Final-turn λ distribution
    # To get the distribution we need to re-run (accumulators only store sums).
    # Instead report the mean final λ and efficiency.
    final_eff  = sum_eff[T_TURNS]  / n
    final_std  = sqrt(max(0.0, sum_eff2[T_TURNS]/n - final_eff^2))
    final_λc   = sum_λ_c[T_TURNS]  / n
    final_λf   = sum_λ_f[T_TURNS]  / n

    @printf("  Final state (t=%d):\n", T_TURNS)
    @printf("    Mean λ_consumers = %.3f\n", final_λc)
    @printf("    Mean λ_firms     = %.3f\n", final_λf)
    @printf("    Efficiency:  mean=%.4f  std=%.4f\n", final_eff, final_std)
    @printf("    vs static Std:    %+.4f\n", final_eff - mean(static_eff_std))
    @printf("    vs static BiGeo:  %+.4f\n", final_eff - mean(static_eff_bgip))
    println()

    (sum_λ_c ./ n, sum_λ_f ./ n, sum_eff ./ n, n)
end

# ── Run all sizes ─────────────────────────────────────────────────────────────

println("ADAPTIVE ZIT — Learning Rule: play whichever strategy last performed better")
println("Performance signal: pro-rata realized trade quantity")
println("Cold start: all agents begin :std; NaN memory forces BiGeo exploration on turn 2")
println()

# Report at these turns
RT = [1, 2, 3, 5, 10, 25, 50, 100, 200, 300]

# Small / mixed markets
let
    # Custom generation loop (heterogeneous sizes, same as zi_efficiency_analysis.jl)
    SEED_BASE = 20_000
    sum_λ_c  = zeros(T_TURNS); sum_λ_f  = zeros(T_TURNS)
    sum_eff  = zeros(T_TURNS); sum_eff2 = zeros(T_TURNS)
    n = 0
    static_std  = Float64[]; static_bgip = Float64[]

    for seed in 1:1_000
        rng   = MersenneTwister(seed + SEED_BASE)
        n_c   = rand(rng, 1:8); n_f = rand(rng, 1:8); max_u = rand(rng, 1:6)
        m, r  = generate_good_market(rng;
            good=1, n_consumers=n_c, n_firms=n_f, max_units=max_u, Q=100)
        r.cleared || continue
        p_star = r.price
        s_wge  = surplus_f64(m, aggregate_demand(m, p_star))
        s_wge > 0 || continue

        consumers = [AdaptConsumer(d, p_star) for d in m.consumers]
        firms     = [AdaptFirm(f, p_star)     for f in m.firms]

        λ_c, λ_f, effs = simulate_market(consumers, firms, m, s_wge; seed=seed)
        for t in 1:T_TURNS
            sum_λ_c[t] += λ_c[t]; sum_λ_f[t] += λ_f[t]
            sum_eff[t]  += effs[t]; sum_eff2[t] += effs[t]^2
        end
        n += 1

        rng2 = MersenneTwister(seed + SEED_BASE + 99_000)
        push!(static_std, mean(begin
            dc = sum(draw_std_consumer(consumers[i], rng2) for i in 1:n_c; init=0)
            ds = sum(draw_std_firm(firms[j], rng2) for j in 1:n_f; init=0)
            surplus_f64(m, min(dc, ds)) / s_wge end for _ in 1:T_TURNS))
        rng3 = MersenneTwister(seed + SEED_BASE + 199_000)
        push!(static_bgip, mean(begin
            dc = sum(draw_bgip_consumer(consumers[i], rng3) for i in 1:n_c; init=0)
            ds = sum(draw_bgip_firm(firms[j], rng3) for j in 1:n_f; init=0)
            surplus_f64(m, min(dc, ds)) / s_wge end for _ in 1:T_TURNS))
    end

    println("=" ^ 66)
    println("Small/mixed markets  (1–8C, 1–8F, 1–6 max_units)  n=$n / 1000")
    @printf("  T_TURNS=%d\n\n", T_TURNS)
    @printf("  Static λ=0 (pure Std):   mean eff = %.4f\n", mean(static_std))
    @printf("  Static λ=1 (pure BiGeo): mean eff = %.4f\n", mean(static_bgip))
    println()
    println("  Turn │  λ_c    λ_f    Eff_mean  Eff_std  vs Std   vs BiGeo")
    println("  ─────┼────────────────────────────────────────────────────────")
    for t in RT
        t > T_TURNS && continue
        λc  = sum_λ_c[t] / n;  λf = sum_λ_f[t] / n
        eff = sum_eff[t]  / n;  sd = sqrt(max(0.0, sum_eff2[t]/n - eff^2))
        @printf("  t=%3d │  %.3f  %.3f  %.4f    %.4f   %+.4f  %+.4f\n",
                t, λc, λf, eff, sd, eff-mean(static_std), eff-mean(static_bgip))
    end
    println()
    @printf("  Final (t=%d): λ_c=%.3f  λ_f=%.3f  eff=%.4f  std=%.4f\n\n",
            T_TURNS, sum_λ_c[T_TURNS]/n, sum_λ_f[T_TURNS]/n,
            sum_eff[T_TURNS]/n, sqrt(max(0.0,sum_eff2[T_TURNS]/n-(sum_eff[T_TURNS]/n)^2)))
end

# 100-agent markets
run_adaptive(nc=50, nf=50, n_markets=300, seed_base=30_000,
             label="100 agents (50C + 50F)  n=300 markets",
             report_turns=RT)

# 500-agent markets
run_adaptive(nc=250, nf=250, n_markets=150, seed_base=31_000,
             label="500 agents (250C + 250F)  n=150 markets",
             report_turns=RT)
