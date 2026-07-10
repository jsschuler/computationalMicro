# Systematic ZI efficiency analysis across random markets.
#
# For each random single-good market with a valid WGE:
#   - ZI consumers get wealth = p* × (their WTP sequence length), so their
#     budget-constrained max demand matches their WGE capacity at p*.
#   - ZI firms get capacity = their WTA sequence length.
#   - Simulation runs T_SIM periods at the fixed WGE price p*.
#   - Efficiency = E[ZI surplus] / WGE surplus, where ZI surplus uses
#     traded = min(demand, supply) per period.
#
# Run with:
#   julia --project=. test/zi_efficiency_analysis.jl

using Random
using Statistics
using Printf
using DiscreteMarket

const N_MARKETS = 1_000
const T_SIM     = 1_000
const SEED_BASE = 20_000

# ── Collection buffers ────────────────────────────────────────────────────────

struct MarketRecord
    efficiency   :: Float64
    n_consumers  :: Int
    n_firms      :: Int
    max_units    :: Int     # draw bound used for generation
    q_wge        :: Int
    surplus_wge  :: Float64
    zi_mean_q    :: Float64
end

function run_analysis()

records = MarketRecord[]
n_no_equil = 0
n_zero_surplus = 0

# ── Main loop ─────────────────────────────────────────────────────────────────

for seed in 1:N_MARKETS
    rng    = MersenneTwister(seed + SEED_BASE)
    n_c    = rand(rng, 1:8)
    n_f    = rand(rng, 1:8)
    max_u  = rand(rng, 1:6)

    m, r = generate_good_market(rng;
        good=1, n_consumers=n_c, n_firms=n_f, max_units=max_u, Q=100)

    r.cleared || (n_no_equil += 1; continue)
    p_star = r.price

    s_wge = Float64(wge_surplus(m, p_star))
    if s_wge <= 0
        n_zero_surplus += 1
        continue
    end

    # ZI agents calibrated to match WGE agent capacities at p*
    zi_consumers = [ZIConsumer(1, p_star * length(d.wtp)) for d in m.consumers]
    zi_firms     = [ZIFirm(1, length(f.wtac))             for f in m.firms]
    zm = ZIMarket(zi_consumers, zi_firms, 1)

    result = zi_simulate(zm, [p_star], T_SIM; seed=seed)
    traded = result.traded[:, 1]

    eff = zi_efficiency(m, p_star, traded)
    isnan(eff) && continue

    push!(records, MarketRecord(
        eff, n_c, n_f, max_u,
        aggregate_demand(m, p_star),
        s_wge,
        mean(traded)
    ))
end

n_analyzed = length(records)
effs = [r.efficiency for r in records]

# ── Summary statistics ────────────────────────────────────────────────────────

println("=" ^ 60)
println("ZI Efficiency Analysis — $N_MARKETS random markets, $T_SIM trials each")
println("=" ^ 60)
println()
println("Markets with no WGE:       $n_no_equil")
println("Markets with zero surplus: $n_zero_surplus")
println("Markets analyzed:          $n_analyzed")
println()
println("── Overall efficiency (ZI surplus / WGE surplus) ──────────")
println("  Mean:    $(round(mean(effs),   digits=4))")
println("  Median:  $(round(median(effs), digits=4))")
println("  Std:     $(round(std(effs),    digits=4))")
println("  Min:     $(round(minimum(effs),digits=4))")
println("  Max:     $(round(maximum(effs),digits=4))")
println()

# ── Histogram ────────────────────────────────────────────────────────────────

println("── Efficiency distribution ─────────────────────────────────")
bins   = 0.0:0.05:1.05
counts = [count(e -> lo <= e < hi, effs) for (lo, hi) in zip(bins, bins[2:end])]
maxc   = max(1, maximum(counts))
barw   = 40
for (i, c) in enumerate(counts)
    lo  = bins[i]
    hi  = bins[i+1]
    bar = "█" ^ round(Int, barw * c / maxc)
    pct = lpad(round(100 * c / n_analyzed, digits=1), 5)
    @printf("  [%.2f–%.2f) │ %-*s %s%% (%d)\n", lo, hi, barw, bar, pct, c)
end
println()

# ── Breakdown by number of consumers ─────────────────────────────────────────

println("── Breakdown by n_consumers ────────────────────────────────")
for group in [(1:2,"1–2"), (3:5,"3–5"), (6:8,"6–8")]
    rng_g, label = group
    sub = [r.efficiency for r in records if r.n_consumers in rng_g]
    isempty(sub) && continue
    @printf("  n_consumers=%s  n=%d  mean=%.4f  median=%.4f  std=%.4f\n",
            label, length(sub), mean(sub), median(sub), std(sub))
end
println()

# ── Breakdown by n_firms ──────────────────────────────────────────────────────

println("── Breakdown by n_firms ────────────────────────────────────")
for group in [(1:2,"1–2"), (3:5,"3–5"), (6:8,"6–8")]
    rng_g, label = group
    sub = [r.efficiency for r in records if r.n_firms in rng_g]
    isempty(sub) && continue
    @printf("  n_firms=%s     n=%d  mean=%.4f  median=%.4f  std=%.4f\n",
            label, length(sub), mean(sub), median(sub), std(sub))
end
println()

# ── Breakdown by max_units ────────────────────────────────────────────────────

println("── Breakdown by max_units (per-agent draw bound) ───────────")
for group in [(1:2,"1–2"), (3:4,"3–4"), (5:6,"5–6")]
    rng_g, label = group
    sub = [r.efficiency for r in records if r.max_units in rng_g]
    isempty(sub) && continue
    @printf("  max_units=%s   n=%d  mean=%.4f  median=%.4f  std=%.4f\n",
            label, length(sub), mean(sub), median(sub), std(sub))
end
println()

# ── Breakdown by WGE quantity ─────────────────────────────────────────────────

println("── Breakdown by WGE equilibrium quantity q* ────────────────")
q_cutoffs = [(1:2,"1–2"), (3:5,"3–5"), (6:10,"6–10"), (11:typemax(Int),"11+")]
for (rng_g, label) in q_cutoffs
    sub = [r.efficiency for r in records if r.q_wge in rng_g]
    isempty(sub) && continue
    @printf("  q*=%s     n=%d (%.1f%%)  mean=%.4f  median=%.4f  std=%.4f\n",
            lpad(label, 3), length(sub), 100*length(sub)/n_analyzed, mean(sub), median(sub), std(sub))
end
println()

# ── High / low efficiency tails ───────────────────────────────────────────────

println("── Tail counts ─────────────────────────────────────────────")
for thresh in [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    n_above = count(e -> e >= thresh, effs)
    @printf("  eff >= %.2f : %d / %d  (%.1f%%)\n",
            thresh, n_above, n_analyzed, 100*n_above/n_analyzed)
end
println()
println("=" ^ 60)

end # run_analysis

run_analysis()
