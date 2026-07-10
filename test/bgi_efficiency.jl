# BiGeometric ZIT efficiency comparison.
#
# Standard ZIT draws quantities from Uniform{0,...,Q_max} (additive).
# BiGeometric ZIT draws from the log-uniform analog:
#   q = floor(exp(u * log(Q_max + 1))) - 1,  u ~ Uniform[0,1]
# giving P(q=k) ∝ 1/k (Zipf/reciprocal on {0,...,Q_max}).
#
# Median = sqrt(Q_max) - 1  vs  Q_max/2 for standard.
# Mean   ≈ Q_max/log(Q_max) vs  Q_max/2 for standard.
#
# This compares standard vs BiGeo ZIT on the same 1000 random markets
# using the same WGE calibration: wealth = p* × len(wtp), capacity = len(wtac).
#
# Run with:
#   julia --project=. test/bgi_efficiency.jl

using Random
using Statistics
using Printf
using DiscreteMarket

const N_MARKETS = 1_000
const T_SIM     = 1_000
const SEED_BASE = 20_000    # same seeds as zi_efficiency_analysis.jl

# ── BiGeometric quantity draw ─────────────────────────────────────────────────
# Log-uniform on {0,...,max_q}: floor(exp(u * log(max_q+1))) - 1

function bgi_draw(max_q::Int, rng::AbstractRNG) :: Int
    max_q <= 0 && return 0
    u = rand(rng)
    floor(Int, exp(u * log(max_q + 2))) - 1
end

# ── Theoretical mean of log-uniform on {0,...,N} ─────────────────────────────
# P(q=k) = log((k+2)/(k+1)) / log(N+2) for k in {0,...,N}
# Inverse CDF: q = floor(exp(u * log(N+2))) - 1, u ~ Uniform[0,1)

bgi_mean_theoretical(max_q::Int) =
    max_q <= 0 ? 0.0 :
    sum(k * log((k+2)/(k+1)) / log(max_q+2) for k in 0:max_q)

# ── Simulate one period under BiGeo ZIT ──────────────────────────────────────

function bgi_simulate(m::GoodMarket, p_star::Price, T::Int; seed::Int=0) :: Vector{Int}
    rng = MersenneTwister(seed)
    traded = Vector{Int}(undef, T)
    for t in 1:T
        demand = sum(bgi_draw(floor(Int, (p_star * length(d.wtp)) / p_star), rng)
                     for d in m.consumers)
        supply = sum(bgi_draw(length(f.wtac), rng) for f in m.firms)
        traded[t] = min(demand, supply)
    end
    traded
end

# ── Surplus helper (Float64, avoids Rational overflow) ────────────────────────

function surplus_f64(m::GoodMarket, q::Int) :: Float64
    q <= 0 && return 0.0
    all_wtp  = sort([Float64(v) for d in m.consumers for v in d.wtp], rev=true)
    all_wtac = sort([Float64(c) for f in m.firms     for c in f.wtac])
    n = min(q, length(all_wtp), length(all_wtac))
    n == 0 && return 0.0
    sum(all_wtp[1:n]) - sum(all_wtac[1:n])
end

# ── Record structs ────────────────────────────────────────────────────────────

struct CompRecord
    q_wge   :: Int
    eff_std :: Float64
    eff_bgi :: Float64
    mean_q_std :: Float64
    mean_q_bgi :: Float64
end

# ── Main ──────────────────────────────────────────────────────────────────────

function run_comparison()

records      = CompRecord[]
n_no_equil   = 0
n_zero_surp  = 0

for seed in 1:N_MARKETS
    rng   = MersenneTwister(seed + SEED_BASE)
    n_c   = rand(rng, 1:8)
    n_f   = rand(rng, 1:8)
    max_u = rand(rng, 1:6)

    m, r = generate_good_market(rng;
        good=1, n_consumers=n_c, n_firms=n_f, max_units=max_u, Q=100)

    r.cleared || (n_no_equil += 1; continue)
    p_star = r.price

    s_wge = surplus_f64(m, aggregate_demand(m, p_star))
    s_wge > 0 || (n_zero_surp += 1; continue)

    q_star = aggregate_demand(m, p_star)

    # ── Standard ZIT ─────────────────────────────────────────────────────────
    zi_consumers = [ZIConsumer(1, p_star * length(d.wtp)) for d in m.consumers]
    zi_firms     = [ZIFirm(1, length(f.wtac))             for f in m.firms]
    zm = ZIMarket(zi_consumers, zi_firms, 1)
    res_std  = zi_simulate(zm, [p_star], T_SIM; seed=seed)
    traded_std = res_std.traded[:, 1]
    eff_std = mean(surplus_f64(m, t) for t in traded_std) / s_wge

    # ── BiGeometric ZIT ──────────────────────────────────────────────────────
    traded_bgi = bgi_simulate(m, p_star, T_SIM; seed=seed)
    eff_bgi = mean(surplus_f64(m, t) for t in traded_bgi) / s_wge

    (isnan(eff_std) || isnan(eff_bgi)) && continue

    push!(records, CompRecord(
        q_star, eff_std, eff_bgi, mean(traded_std), mean(Float64.(traded_bgi))
    ))
end

n = length(records)
effs_std = [r.eff_std for r in records]
effs_bgi = [r.eff_bgi for r in records]

# ── Print mean-demand comparison per agent for a few Q_max values ─────────────

println("=" ^ 64)
println("BiGeometric ZIT — Log-Uniform Quantity Draw")
println("  q = floor(exp(u * log(Q_max+1))) - 1,  u ~ Uniform[0,1]")
println("=" ^ 64)
println()
println("── Theoretical mean demand by Q_max ──────────────────────────")
println("  Q_max │  Standard   BiGeo    Ratio    Median_bgi")
for qm in [1, 2, 5, 10, 25, 50, 100, 250]
    std_mean = qm / 2.0
    bgi_mean = bgi_mean_theoretical(qm)
    med_bgi  = floor(Int, sqrt(qm + 1)) - 1
    @printf("  %5d │  %6.2f    %6.2f   %5.3f    %d\n",
            qm, std_mean, bgi_mean, bgi_mean/std_mean, med_bgi)
end
println()

# ── Side-by-side efficiency summary ──────────────────────────────────────────

println("── Efficiency comparison (n=$n markets, $T_SIM trials each) ────")
println()
println("                Standard ZIT    BiGeo ZIT    Δ (BiGeo − Std)")
@printf("  Mean:         %.4f          %.4f       %+.4f\n",
        mean(effs_std), mean(effs_bgi), mean(effs_bgi) - mean(effs_std))
@printf("  Median:       %.4f          %.4f       %+.4f\n",
        median(effs_std), median(effs_bgi), median(effs_bgi) - median(effs_std))
@printf("  Std:          %.4f          %.4f       %+.4f\n",
        std(effs_std), std(effs_bgi), std(effs_bgi) - std(effs_std))
@printf("  Min:          %.4f         %.4f      %+.4f\n",
        minimum(effs_std), minimum(effs_bgi), minimum(effs_bgi) - minimum(effs_std))
@printf("  Max:          %.4f          %.4f       %+.4f\n",
        maximum(effs_std), maximum(effs_bgi), maximum(effs_bgi) - maximum(effs_std))
println()

println("── Negative efficiency rate ──────────────────────────────────")
n_neg_std = count(e -> e < 0, effs_std)
n_neg_bgi = count(e -> e < 0, effs_bgi)
@printf("  Standard: %d / %d  (%.1f%%)\n", n_neg_std, n, 100*n_neg_std/n)
@printf("  BiGeo:    %d / %d  (%.1f%%)\n", n_neg_bgi, n, 100*n_neg_bgi/n)
println()

println("── Tail counts ───────────────────────────────────────────────")
println("  Threshold │  Standard    BiGeo")
for thresh in [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    n_std = count(e -> e >= thresh, effs_std)
    n_bgi = count(e -> e >= thresh, effs_bgi)
    @printf("  eff≥%.2f  │  %3d/%d (%.1f%%)  %3d/%d (%.1f%%)\n",
            thresh,
            n_std, n, 100*n_std/n,
            n_bgi, n, 100*n_bgi/n)
end
println()

# ── Breakdown by q* ──────────────────────────────────────────────────────────

println("── Breakdown by WGE quantity q* ──────────────────────────────")
println("  q* range  │  n    Std mean   BiGeo mean   Δ")
for (rng_q, label) in [(1:2,"1–2"), (3:5,"3–5"), (6:10,"6–10"), (11:typemax(Int),"11+")]
    sub = [r for r in records if r.q_wge in rng_q]
    isempty(sub) && continue
    ms = mean(r.eff_std for r in sub)
    mb = mean(r.eff_bgi for r in sub)
    @printf("  q*=%-5s   │  %-4d  %.4f     %.4f      %+.4f\n",
            label, length(sub), ms, mb, mb - ms)
end
println()

# ── Scatter summary: how often does BiGeo beat standard? ─────────────────────

n_bgi_wins  = count(r -> r.eff_bgi > r.eff_std, records)
n_std_wins  = count(r -> r.eff_bgi < r.eff_std, records)
n_tie       = count(r -> r.eff_bgi ≈ r.eff_std, records)
@printf("── Head-to-head (same market, same seed) ─────────────────────\n")
@printf("  BiGeo > Std:  %d / %d  (%.1f%%)\n", n_bgi_wins, n, 100*n_bgi_wins/n)
@printf("  Std > BiGeo:  %d / %d  (%.1f%%)\n", n_std_wins, n, 100*n_std_wins/n)
@printf("  Tie:          %d / %d  (%.1f%%)\n", n_tie,      n, 100*n_tie/n)
mean_improvement = mean(r.eff_bgi - r.eff_std for r in records)
@printf("  Mean BiGeo−Std improvement: %+.4f\n", mean_improvement)
println()

# ── Large-market projection ───────────────────────────────────────────────────

println("── Large-market scaling prediction ───────────────────────────")
println("  With max k units/agent:")
println("    Standard E[q_i] = k/2")
println("    BiGeo E[q_i]    = sum_{j=0}^{k} j*log((j+2)/(j+1)) / log(k+2)")
println()
println("  k    │  Std E[q]  BiGeo E[q]  Ratio   Predicted BiGeo eff. (large N)")
for k in [3, 5, 10, 25, 50]
    std_mean = k / 2.0
    bgi_mean = bgi_mean_theoretical(k)
    ratio    = bgi_mean / std_mean
    @printf("  k=%-3d │  %6.2f     %6.2f     %.3f   %.1f%%\n",
            k, std_mean, bgi_mean, ratio, 100*ratio)
end
println()
println("  (Standard ZIT → 100% as N→∞; BiGeo ZIT → const < 100%.)")
println()

end # run_comparison

run_comparison()
