# Adaptive ZIT with individual SURPLUS as the performance signal.
#
# Identical learning rule and structure as adaptive_zit.jl, but agents observe
# their own consumer/producer surplus from the units they actually trade —
# not the raw number of units realized.
#
# Individual surplus definitions:
#   Consumer i assigned k units:
#     surplus = sum(active_wtp[1:min(k, n_active)]) - min(k, n_active) * p*
#   Firm j assigned k units:
#     surplus = min(k, n_active) * p* - sum(active_wtac[1:min(k, n_active)])
#
# Both are capped at n_active: units with wtp < p* or wtac > p* contribute 0.
# cum_surplus[k] is precomputed so the lookup is O(1).
#
# Prediction: BiGeo price-space selects high-wtp inframarginal units, giving
# higher per-unit surplus.  Surplus signal should let thin-market consumers
# identify BiGeo as superior, unlike the quantity signal.
#
# Run with:
#   julia --project=. test/adaptive_zit_surplus.jl

using Random
using Statistics
using Printf
using DiscreteMarket

const T_TURNS   = 300
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

# ── Precomputed per-agent data ─────────────────────────────────────────────────
# cum_surplus[k] = individual surplus from trading k units (k=0 not stored; use 0.0)
# cum_surplus[k] for k ∈ 1:n_active is precomputed; values above n_active repeat [n_active].

struct AdaptConsumer
    n_total     :: Int
    log_p       :: Float64
    log_hi      :: Float64
    active      :: Vector{Float64}   # wtp ≥ p*, sorted descending
    cum_surplus :: Vector{Float64}   # cum_surplus[k] = sum(active[1:k]) - k*p*, length = n_active
end

struct AdaptFirm
    n_total     :: Int
    log_lo      :: Float64
    log_p       :: Float64
    active      :: Vector{Float64}   # wtac ≤ p*, sorted ascending
    cum_surplus :: Vector{Float64}   # cum_surplus[k] = k*p* - sum(active[1:k]), length = n_active
end

function AdaptConsumer(d::ConsumerDemand, p_star::Price)
    lp  = log(Float64(p_star))
    act = sort([Float64(v) for v in d.wtp if v >= p_star], rev=true)
    lhi = isempty(act) ? lp : log(act[1])
    pf  = Float64(p_star)
    cs  = isempty(act) ? Float64[] :
          [sum(act[1:k]) - k * pf for k in 1:length(act)]
    AdaptConsumer(length(d.wtp), lp, lhi, act, cs)
end

function AdaptFirm(f::FirmSupply, p_star::Price)
    lp  = log(Float64(p_star))
    act = sort([Float64(c) for c in f.wtac if c <= p_star])
    llo = isempty(act) ? lp : log(act[1])
    pf  = Float64(p_star)
    cs  = isempty(act) ? Float64[] :
          [k * pf - sum(act[1:k]) for k in 1:length(act)]
    AdaptFirm(length(f.wtac), llo, lp, act, cs)
end

# ── Individual surplus lookup ─────────────────────────────────────────────────

@inline function ind_surplus_c(c::AdaptConsumer, k::Int) :: Float64
    k <= 0 && return 0.0
    isempty(c.cum_surplus) && return 0.0
    c.cum_surplus[min(k, length(c.cum_surplus))]
end

@inline function ind_surplus_f(f::AdaptFirm, k::Int) :: Float64
    k <= 0 && return 0.0
    isempty(f.cum_surplus) && return 0.0
    f.cum_surplus[min(k, length(f.cum_surplus))]
end

# ── Single-agent draws ────────────────────────────────────────────────────────

@inline draw_std_consumer(c::AdaptConsumer, rng::AbstractRNG) :: Int =
    rand(rng, 0:c.n_total)

@inline function draw_bgip_consumer(c::AdaptConsumer, rng::AbstractRNG) :: Int
    isempty(c.active) && return 0
    c.log_hi ≈ c.log_p && return length(c.active)
    pv = exp(c.log_p + rand(rng) * (c.log_hi - c.log_p))
    searchsortedlast(c.active, pv, rev=true)
end

@inline draw_std_firm(f::AdaptFirm, rng::AbstractRNG) :: Int =
    rand(rng, 0:f.n_total)

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

@inline function next_strategy(perf_std::Float64, perf_bgip::Float64) :: Bool
    ps = isnan(perf_std)  ? +Inf : perf_std
    pb = isnan(perf_bgip) ? +Inf : perf_bgip
    pb > ps
end

# ── Simulate one market for T_TURNS turns ────────────────────────────────────

function simulate_market(consumers::Vector{AdaptConsumer},
                          firms::Vector{AdaptFirm},
                          m::GoodMarket, s_wge::Float64;
                          seed::Int=0) :: Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}

    nc, nf = length(consumers), length(firms)

    use_bgip_c   = falses(nc)
    perf_std_c   = fill(NaN, nc)
    perf_bgip_c  = fill(NaN, nc)

    use_bgip_f   = falses(nf)
    perf_std_f   = fill(NaN, nf)
    perf_bgip_f  = fill(NaN, nf)

    λ_c  = Vector{Float64}(undef, T_TURNS)
    λ_f  = Vector{Float64}(undef, T_TURNS)
    effs = Vector{Float64}(undef, T_TURNS)

    rng     = MersenneTwister(seed)
    draws_c = Vector{Int}(undef, nc)
    draws_f = Vector{Int}(undef, nf)

    for t in 1:T_TURNS
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

        # Performance update: individual surplus from pro-rata realized units
        for i in 1:nc
            r = realized_c(draws_c[i], D, S)
            sig = ind_surplus_c(consumers[i], r)
            if use_bgip_c[i]
                perf_bgip_c[i] = sig
            else
                perf_std_c[i]  = sig
            end
        end
        for j in 1:nf
            r = realized_f(draws_f[j], D, S)
            sig = ind_surplus_f(firms[j], r)
            if use_bgip_f[j]
                perf_bgip_f[j] = sig
            else
                perf_std_f[j]  = sig
            end
        end

        for i in 1:nc
            use_bgip_c[i] = next_strategy(perf_std_c[i], perf_bgip_c[i])
        end
        for j in 1:nf
            use_bgip_f[j] = next_strategy(perf_std_f[j], perf_bgip_f[j])
        end
    end

    λ_c, λ_f, effs
end

# ── Run across many markets ───────────────────────────────────────────────────

function run_adaptive(; nc, nf, n_markets, seed_base, label, max_units=MAX_UNITS,
                       report_turns=[1,2,3,5,10,25,50,100,200,300])

    sum_λ_c  = zeros(T_TURNS); sum_λ_f  = zeros(T_TURNS)
    sum_eff  = zeros(T_TURNS); sum_eff2 = zeros(T_TURNS)
    n_analyzed = 0
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
            sum_λ_c[t]  += λ_c[t]; sum_λ_f[t]  += λ_f[t]
            sum_eff[t]  += effs[t]; sum_eff2[t] += effs[t]^2
        end
        n_analyzed += 1

        rng2 = MersenneTwister(seed + seed_base + 99_000)
        push!(static_eff_std, mean(begin
            dc = sum(draw_std_consumer(consumers[i], rng2) for i in 1:nc; init=0)
            ds = sum(draw_std_firm(firms[j], rng2) for j in 1:nf; init=0)
            surplus_f64(m, min(dc, ds)) / s_wge
        end for _ in 1:T_TURNS))

        rng3 = MersenneTwister(seed + seed_base + 199_000)
        push!(static_eff_bgip, mean(begin
            dc = sum(draw_bgip_consumer(consumers[i], rng3) for i in 1:nc; init=0)
            ds = sum(draw_bgip_firm(firms[j], rng3) for j in 1:nf; init=0)
            surplus_f64(m, min(dc, ds)) / s_wge
        end for _ in 1:T_TURNS))
    end

    n = n_analyzed
    println("=" ^ 66)
    println(label)
    @printf("  n_analyzed=%d / %d  |  T_TURNS=%d\n\n", n, n_markets, T_TURNS)
    @printf("  Static λ=0 (pure Std):   mean eff = %.4f\n", mean(static_eff_std))
    @printf("  Static λ=1 (pure BiGeo): mean eff = %.4f\n", mean(static_eff_bgip))
    println()

    println("  Turn │  λ_c    λ_f    Eff_mean  Eff_std  vs Std   vs BiGeo")
    println("  ─────┼────────────────────────────────────────────────────────")
    for t in report_turns
        t > T_TURNS && continue
        λc  = sum_λ_c[t] / n;  λf  = sum_λ_f[t] / n
        eff = sum_eff[t]  / n;  var = sum_eff2[t]/n - eff^2
        sd  = var > 0 ? sqrt(var) : 0.0
        @printf("  t=%3d │  %.3f  %.3f  %.4f    %.4f   %+.4f  %+.4f\n",
                t, λc, λf, eff, sd,
                eff - mean(static_eff_std), eff - mean(static_eff_bgip))
    end
    println()

    final_eff = sum_eff[T_TURNS]  / n
    final_std = sqrt(max(0.0, sum_eff2[T_TURNS]/n - final_eff^2))
    @printf("  Final state (t=%d):\n", T_TURNS)
    @printf("    Mean λ_consumers = %.3f\n", sum_λ_c[T_TURNS] / n)
    @printf("    Mean λ_firms     = %.3f\n", sum_λ_f[T_TURNS] / n)
    @printf("    Efficiency:  mean=%.4f  std=%.4f\n", final_eff, final_std)
    @printf("    vs static Std:    %+.4f\n", final_eff - mean(static_eff_std))
    @printf("    vs static BiGeo:  %+.4f\n", final_eff - mean(static_eff_bgip))
    println()

    (sum_λ_c ./ n, sum_λ_f ./ n, sum_eff ./ n, n)
end

# ── Main ──────────────────────────────────────────────────────────────────────

println("ADAPTIVE ZIT — Learning Rule: play whichever strategy last performed better")
println("Performance signal: individual surplus from pro-rata realized units")
println("  Consumer surplus(k) = sum(active_wtp[1:k]) - k*p*  (capped at n_active)")
println("  Firm     surplus(k) = k*p* - sum(active_wtac[1:k]) (capped at n_active)")
println()

RT = [1, 2, 3, 5, 10, 25, 50, 100, 200, 300]

# Small / mixed markets
let
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
            sum_eff[T_TURNS]/n,
            sqrt(max(0.0, sum_eff2[T_TURNS]/n - (sum_eff[T_TURNS]/n)^2)))
end

# 100-agent markets
run_adaptive(nc=50, nf=50, n_markets=300, seed_base=30_000,
             label="100 agents (50C + 50F)  n=300 markets",
             report_turns=RT)

# 500-agent markets
run_adaptive(nc=250, nf=250, n_markets=150, seed_base=31_000,
             label="500 agents (250C + 250F)  n=150 markets",
             report_turns=RT)
