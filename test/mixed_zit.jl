# Mixture model: additive (standard) ZIT + BiGeometric price-space ZIT.
#
# Mixing parameter λ ∈ [0,1].  Each agent, each period, independently:
#   with prob (1-λ): draw via standard ZIT  (q ~ Uniform{0,...,Q_max})
#   with prob    λ:  draw via BiGeo price   (p_v ~ LogUniform[p*, max_WTP])
#
# λ=0 → pure standard ZIT; λ=1 → pure BiGeo price-space ZIT.
# Sweeping λ tests whether efficiency interpolates smoothly between the two.
#
# Run with:
#   julia --project=. test/mixed_zit.jl

using Random
using Statistics
using Printf
using DiscreteMarket

const T_SIM = 500   # trials per market (reduced; we're sweeping λ)

# ── Surplus ───────────────────────────────────────────────────────────────────

function surplus_f64(m::GoodMarket, q::Int) :: Float64
    q <= 0 && return 0.0
    all_wtp  = sort([Float64(v) for d in m.consumers for v in d.wtp], rev=true)
    all_wtac = sort([Float64(c) for f in m.firms     for c in f.wtac])
    n = min(q, length(all_wtp), length(all_wtac))
    n == 0 && return 0.0
    sum(all_wtp[1:n]) - sum(all_wtac[1:n])
end

# ── Precomputed per-agent data for the mixed draw ─────────────────────────────

struct MixedConsumer
    n_total    :: Int           # standard ZIT max_q = floor(wealth/p*) = len(wtp)
    log_p      :: Float64       # log(p*)
    log_max    :: Float64       # log(max active WTP); = log_p if no spread
    active_wtp :: Vector{Float64}  # sorted descending, WTP ≥ p*
end

struct MixedFirm
    n_total    :: Int           # standard ZIT capacity = len(wtac)
    log_min    :: Float64       # log(min active WTAC); = log_p if no spread
    log_p      :: Float64       # log(p*)
    active_wtac :: Vector{Float64} # sorted ascending, WTAC ≤ p*
end

function MixedConsumer(d::ConsumerDemand, p_star::Price)
    lp     = log(Float64(p_star))
    active = sort([Float64(v) for v in d.wtp if v >= p_star], rev=true)
    lhi    = isempty(active) ? lp : log(active[1])
    MixedConsumer(length(d.wtp), lp, lhi, active)
end

function MixedFirm(f::FirmSupply, p_star::Price)
    lp     = log(Float64(p_star))
    active = sort([Float64(c) for c in f.wtac if c <= p_star])
    llo    = isempty(active) ? lp : log(active[1])
    MixedFirm(length(f.wtac), llo, lp, active)
end

# ── Mixed draw ────────────────────────────────────────────────────────────────

@inline function mixed_demand(c::MixedConsumer, λ::Float64, rng::AbstractRNG) :: Int
    if rand(rng) < λ
        # BiGeo price: p_v ~ LogUniform[p*, max_wtp]
        isempty(c.active_wtp) && return 0
        c.log_max ≈ c.log_p && return length(c.active_wtp)
        pv = exp(c.log_p + rand(rng) * (c.log_max - c.log_p))
        return searchsortedlast(c.active_wtp, pv, rev=true)
    else
        # Standard ZIT: q ~ Uniform{0,...,n_total}
        return rand(rng, 0:c.n_total)
    end
end

@inline function mixed_supply(f::MixedFirm, λ::Float64, rng::AbstractRNG) :: Int
    if rand(rng) < λ
        # BiGeo price: p_v ~ LogUniform[min_wtac, p*]
        isempty(f.active_wtac) && return 0
        f.log_p ≈ f.log_min && return length(f.active_wtac)
        pv = exp(f.log_min + rand(rng) * (f.log_p - f.log_min))
        return searchsortedlast(f.active_wtac, pv)
    else
        # Standard ZIT: q ~ Uniform{0,...,n_total}
        return rand(rng, 0:f.n_total)
    end
end

# ── Simulate one market across all λ values in a single pass ─────────────────
# Pre-draws T_SIM × n_agents × 2 uniforms; applies each λ post-hoc.
# This makes λ values comparable on the same random seeds.

function simulate_λ_sweep(consumers::Vector{MixedConsumer},
                           firms::Vector{MixedFirm},
                           λ_vals::Vector{Float64},
                           s_wge::Float64,
                           m::GoodMarket;
                           seed::Int=0) :: Vector{Float64}
    nc, nf = length(consumers), length(firms)
    # Pre-draw: coin flip + draw value, per agent per period
    rng      = MersenneTwister(seed)
    coins_c  = rand(rng, nc, T_SIM)   # coin flip for each consumer×period
    draws_c  = rand(rng, nc, T_SIM)   # draw value for each consumer×period
    coins_f  = rand(rng, nf, T_SIM)
    draws_f  = rand(rng, nf, T_SIM)

    effs = Vector{Float64}(undef, length(λ_vals))

    for (li, λ) in enumerate(λ_vals)
        total_surplus = 0.0
        for t in 1:T_SIM
            demand = 0
            for i in 1:nc
                c = consumers[i]
                if coins_c[i,t] < λ
                    # BiGeo price
                    if !isempty(c.active_wtp) && !(c.log_max ≈ c.log_p)
                        pv = exp(c.log_p + draws_c[i,t] * (c.log_max - c.log_p))
                        demand += searchsortedlast(c.active_wtp, pv, rev=true)
                    else
                        demand += length(c.active_wtp)
                    end
                else
                    # Standard ZIT
                    demand += floor(Int, draws_c[i,t] * (c.n_total + 1))
                end
            end
            supply = 0
            for i in 1:nf
                f = firms[i]
                if coins_f[i,t] < λ
                    # BiGeo price
                    if !isempty(f.active_wtac) && !(f.log_p ≈ f.log_min)
                        pv = exp(f.log_min + draws_f[i,t] * (f.log_p - f.log_min))
                        supply += searchsortedlast(f.active_wtac, pv)
                    else
                        supply += length(f.active_wtac)
                    end
                else
                    # Standard ZIT
                    supply += floor(Int, draws_f[i,t] * (f.n_total + 1))
                end
            end
            total_surplus += surplus_f64(m, min(demand, supply))
        end
        effs[li] = (total_surplus / T_SIM) / s_wge
    end

    effs
end

# ── Run across markets ────────────────────────────────────────────────────────

function run_sweep(; nc, nf, n_markets, seed_base, λ_vals, max_units=5, Q=100,
                   label="")
    all_effs = [Float64[] for _ in λ_vals]
    n_skip   = 0

    for seed in 1:n_markets
        rng = MersenneTwister(seed + seed_base)
        m, r = generate_good_market(rng;
            good=1, n_consumers=nc, n_firms=nf,
            max_units=max_units, Q=Q)
        r.cleared || (n_skip += 1; continue)
        p_star = r.price
        s_wge  = surplus_f64(m, aggregate_demand(m, p_star))
        s_wge > 0 || (n_skip += 1; continue)

        consumers = [MixedConsumer(d, p_star) for d in m.consumers]
        firms     = [MixedFirm(f, p_star)     for f in m.firms]

        effs = simulate_λ_sweep(consumers, firms, λ_vals, s_wge, m; seed=seed)
        for (li, e) in enumerate(effs)
            isnan(e) || push!(all_effs[li], e)
        end
    end

    n_analyzed = length(all_effs[1])
    println("=" ^ 62)
    println(label == "" ? "$nc consumers + $nf firms  ($(nc+nf) agents)" : label)
    @printf("  n_markets=%d  T_SIM=%d  skipped=%d  analyzed=%d\n\n",
            n_markets, T_SIM, n_skip, n_analyzed)

    println("  λ     │  Mean    Median   Std     Min     >0.90   >0.95")
    println("  ──────┼────────────────────────────────────────────────")
    for (li, λ) in enumerate(λ_vals)
        es = all_effs[li]
        isempty(es) && continue
        n90 = count(e -> e >= 0.90, es)
        n95 = count(e -> e >= 0.95, es)
        @printf("  λ=%.2f │  %.4f  %.4f  %.4f  %7.4f  %5.1f%%  %5.1f%%\n",
                λ, mean(es), median(es), std(es), minimum(es),
                100*n90/length(es), 100*n95/length(es))
    end
    println()

    # Head-to-head: each λ vs λ=0 (standard ZIT baseline)
    println("  Head-to-head vs λ=0 (standard ZIT):")
    base = all_effs[1]
    for (li, λ) in enumerate(λ_vals)
        li == 1 && continue
        es   = all_effs[li]
        n    = min(length(base), length(es))
        wins = count(i -> es[i] > base[i], 1:n)
        @printf("  λ=%.2f: wins %d/%d (%.1f%%)  mean Δ=%+.4f\n",
                λ, wins, n, 100*wins/n, mean(es[1:n]) - mean(base[1:n]))
    end
    println()

    all_effs
end

# ── Small / mixed markets ─────────────────────────────────────────────────────

λ_fine = collect(0.0:0.1:1.0)

println("SMALL / MIXED MARKETS  (1–8 consumers, 1–8 firms, 1–6 max_units)")
println()

# Pre-run: use same market-generation logic as zi_efficiency_analysis.jl
small_effs = let
    all_effs = [Float64[] for _ in λ_fine]
    n_skip   = 0
    SEED_BASE = 20_000

    for seed in 1:1_000
        rng   = MersenneTwister(seed + SEED_BASE)
        n_c   = rand(rng, 1:8)
        n_f   = rand(rng, 1:8)
        max_u = rand(rng, 1:6)
        m, r  = generate_good_market(rng;
            good=1, n_consumers=n_c, n_firms=n_f, max_units=max_u, Q=100)
        r.cleared || (n_skip += 1; continue)
        p_star = r.price
        s_wge  = surplus_f64(m, aggregate_demand(m, p_star))
        s_wge > 0 || (n_skip += 1; continue)

        consumers = [MixedConsumer(d, p_star) for d in m.consumers]
        firms     = [MixedFirm(f, p_star)     for f in m.firms]
        effs = simulate_λ_sweep(consumers, firms, λ_fine, s_wge, m; seed=seed)
        for (li, e) in enumerate(effs)
            isnan(e) || push!(all_effs[li], e)
        end
    end

    n_analyzed = length(all_effs[1])
    println("=" ^ 62)
    println("Small/mixed markets (1–8C, 1–8F, 1–6 max_units)  n=1000")
    @printf("  T_SIM=%d  skipped=%d  analyzed=%d\n\n", T_SIM, n_skip, n_analyzed)

    println("  λ     │  Mean    Median   Std     Min     >0.50   >0.90")
    println("  ──────┼────────────────────────────────────────────────")
    for (li, λ) in enumerate(λ_fine)
        es = all_effs[li]
        isempty(es) && continue
        n50 = count(e -> e >= 0.50, es)
        n90 = count(e -> e >= 0.90, es)
        @printf("  λ=%.1f │  %.4f  %.4f  %.4f  %7.4f  %5.1f%%  %5.1f%%\n",
                λ, mean(es), median(es), std(es), minimum(es),
                100*n50/length(es), 100*n90/length(es))
    end
    println()

    base = all_effs[1]
    println("  Head-to-head vs λ=0 (standard ZIT):")
    for (li, λ) in enumerate(λ_fine)
        li == 1 && continue
        es   = all_effs[li]
        n    = min(length(base), length(es))
        wins = count(i -> es[i] > base[i], 1:n)
        @printf("  λ=%.1f: wins %d/%d (%.1f%%)  mean Δ=%+.4f\n",
                λ, wins, n, 100*wins/n, mean(es[1:n]) - mean(base[1:n]))
    end
    println()

    all_effs
end

# ── 100-agent markets ─────────────────────────────────────────────────────────

println("100-AGENT MARKETS  (50C + 50F, max 5 units/agent)")
println()
λ_coarse = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
run_sweep(nc=50, nf=50, n_markets=300, seed_base=30_000, λ_vals=λ_coarse,
          label="100 agents (50C + 50F)  n=300 markets")

# ── 500-agent markets ─────────────────────────────────────────────────────────

println("500-AGENT MARKETS  (250C + 250F, max 5 units/agent)")
println()
run_sweep(nc=250, nf=250, n_markets=200, seed_base=31_000, λ_vals=λ_coarse,
          label="500 agents (250C + 250F)  n=200 markets")
