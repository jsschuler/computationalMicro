# ZI vs WGE comparison in large markets.
# Runs two sizes: 100 agents (50C+50F) and 500 agents (250C+250F).
#
# Run with:
#   julia --project=. test/zi_large_market.jl

using Random
using Statistics
using Printf
using DiscreteMarket

const N_MARKETS = 500
const T_SIM     = 1_000
const MAX_UNITS = 5

# total_surplus uses Rational{Int64} which overflows at large market sizes.
# This float version is safe for analysis; exact rational is only needed for tests.
function surplus_f64(m::GoodMarket, q::Int) :: Float64
    q <= 0 && return 0.0
    all_wtp  = sort([Float64(v) for d in m.consumers for v in d.wtp],  rev=true)
    all_wtac = sort([Float64(c) for f in m.firms     for c in f.wtac])
    n = min(q, length(all_wtp), length(all_wtac))
    n == 0 && return 0.0
    sum(all_wtp[1:n]) - sum(all_wtac[1:n])
end

struct LargeRecord
    q_wge        :: Int
    surplus_wge  :: Float64
    zi_mean_q    :: Float64
    zi_mean_surp :: Float64
    efficiency   :: Float64
end

function run_size(n_consumers, n_firms; seed_base)
    records    = LargeRecord[]
    n_no_equil = 0

    for seed in 1:N_MARKETS
        rng = MersenneTwister(seed + seed_base)
        m, r = generate_good_market(rng;
            good=1, n_consumers=n_consumers, n_firms=n_firms,
            max_units=MAX_UNITS, Q=100)

        r.cleared || (n_no_equil += 1; continue)
        p_star = r.price
        s_wge  = surplus_f64(m, aggregate_demand(m, p_star))
        s_wge > 0 || continue

        zi_consumers = [ZIConsumer(1, p_star * length(d.wtp)) for d in m.consumers]
        zi_firms     = [ZIFirm(1, length(f.wtac))             for f in m.firms]
        zm = ZIMarket(zi_consumers, zi_firms, 1)

        result = zi_simulate(zm, [p_star], T_SIM; seed=seed)
        traded = result.traded[:, 1]

        zi_mean_surp = mean(surplus_f64(m, t) for t in traded)
        eff = zi_mean_surp / s_wge
        isnan(eff) && continue

        push!(records, LargeRecord(
            aggregate_demand(m, p_star), s_wge,
            mean(traded), Float64(zi_mean_surp), eff))
    end

    records, n_no_equil
end

function print_size(n_consumers, n_firms, records, n_no_equil)
    n     = length(records)
    effs  = [r.efficiency   for r in records]
    q_wge = [r.q_wge        for r in records]
    q_zi  = [r.zi_mean_q    for r in records]
    s_wge = [r.surplus_wge  for r in records]
    s_zi  = [r.zi_mean_surp for r in records]

    println("=" ^ 62)
    println("$n_consumers consumers + $n_firms firms  " *
            "($(n_consumers+n_firms) agents total)")
    println("  max $MAX_UNITS units/agent, $N_MARKETS seeds, $T_SIM ZI trials each")
    println("=" ^ 62)
    @printf("Markets analyzed:  %d / %d  (%d had no WGE)\n\n", n, N_MARKETS, n_no_equil)

    println("── WGE outcomes ────────────────────────────────────────────")
    @printf("  q*:  mean=%.1f  median=%.1f  [%.0f, %.0f]\n",
            mean(q_wge), median(q_wge), minimum(q_wge), maximum(q_wge))
    @printf("  S*:  mean=%.1f  median=%.1f  [%.1f, %.1f]\n\n",
            mean(s_wge), median(s_wge), minimum(s_wge), maximum(s_wge))

    println("── ZI outcomes (at fixed p*) ───────────────────────────────")
    @printf("  qty: mean=%.1f  median=%.1f  [%.1f, %.1f]\n",
            mean(q_zi), median(q_zi), minimum(q_zi), maximum(q_zi))
    @printf("  S:   mean=%.1f  median=%.1f  [%.1f, %.1f]\n\n",
            mean(s_zi), median(s_zi), minimum(s_zi), maximum(s_zi))

    @printf("── WGE vs ZI gap\n")
    @printf("  Mean quantity gap  (q* − ZI):  %.2f units\n",  mean(q_wge .- q_zi))
    @printf("  Mean surplus ratio (ZI / WGE): %.4f\n\n", mean(s_zi ./ s_wge))

    println("── ZI efficiency ───────────────────────────────────────────")
    @printf("  Mean=%.4f  Median=%.4f  Std=%.4f  Min=%.4f  Max=%.4f\n\n",
            mean(effs), median(effs), std(effs), minimum(effs), maximum(effs))

    println("── Tail counts ─────────────────────────────────────────────")
    for thresh in [0.90, 0.95, 0.98, 0.99]
        n_above = count(e -> e >= thresh, effs)
        @printf("  eff >= %.2f : %d / %d  (%.1f%%)\n",
                thresh, n_above, n, 100*n_above/n)
    end
    println()
end

for (nc, nf, sb) in [(50, 50, 30_000), (250, 250, 31_000)]
    recs, no_eq = run_size(nc, nf; seed_base=sb)
    print_size(nc, nf, recs, no_eq)
end
