# BiGeometric jump index test.
#
# Two predictors of WGE non-existence are compared against J (additive):
#
#   BJ(m)  = max_v  J(v) / D(v+)    where D(v+) = D(v) - J(v) > 0
#            "relative jump from the perspective of post-jump demand"
#            Excludes the global maximum WTP where D(v+) = 0 (trivially 1).
#
#   BJ2(m) = max_v  J(v) / z(v)     where z(v) = D(v) - S_max(v) > 0
#            "how much the jump overshoots the supply gap"
#            Only defined at prices with excess demand.
#            BJ2 > 1 is the direct BiGeometric bridge-failure condition.
#
# Also tests the exact bridge-failure condition (see below) and checks
# whether BJ2 > 1 ↔ bridge_failure ↔ WGE non-existence.
#
# Run with:
#   julia --project=. test/bj_index.jl

using Random
using Statistics
using Printf
using DiscreteMarket

# ── Additive jump index ───────────────────────────────────────────────────────

function jump_index(m::GoodMarket) :: Int
    counts = Dict{Price, Int}()
    for d in m.consumers, v in d.wtp
        counts[v] = get(counts, v, 0) + 1
    end
    isempty(counts) ? 0 : maximum(values(counts))
end

# ── BJ: relative jump vs post-jump demand ─────────────────────────────────────
# BJ(v) = J(v) / D(v+) where D(v+) = D(v) - J(v).
# Skips v_max where D(v+) = 0.

function bj_index(m::GoodMarket) :: Float64
    all_wtp = vcat([d.wtp for d in m.consumers]...)
    isempty(all_wtp) && return 0.0
    bj = 0.0
    for v in unique(all_wtp)
        j_v   = count(==(v), all_wtp)
        d_v   = aggregate_demand(m, v)
        d_after = d_v - j_v
        d_after > 0 || continue          # skip max WTP (D(v+) = 0)
        bj = max(bj, j_v / d_after)
    end
    bj
end

# ── BJ2: jump vs supply gap — the direct bridge-failure ratio ─────────────────
# BJ2(v) = J(v) / z(v)  where z(v) = D(v) - S_max(v) > 0.
# BJ2 > 1 at some v ↔ jump overshoots supply → bridge failure.

function bj2_index(m::GoodMarket) :: Float64
    all_wtp = vcat([d.wtp for d in m.consumers]...)
    isempty(all_wtp) && return 0.0
    bj2 = 0.0
    for v in unique(all_wtp)
        j_v = count(==(v), all_wtp)
        j_v > 1 || continue              # J=1 can never overshoot (z ≥ 1)
        z_v = excess_demand(m, v)
        z_v > 0 || continue              # only at excess-demand prices
        bj2 = max(bj2, j_v / z_v)
    end
    bj2
end

# ── Exact bridge-failure test ─────────────────────────────────────────────────
# For each WTP value v (iterated highest-to-lowest):
#   excess_at_v    = D(v) - S_max(v)          > 0 → excess demand at v
#   excess_above_v = D(v)-J(v) - S_max(v)     < 0 → excess supply just above v
# When both hold, the price line is partitioned: excess demand at all p ≤ v and
# excess supply at all p > v (supply non-decreasing, demand non-increasing).
# No clearing exists anywhere.  This condition is equivalent to BJ2 > 1.

function bridge_failure(m::GoodMarket) :: Bool
    all_wtp_raw = vcat([d.wtp for d in m.consumers]...)
    isempty(all_wtp_raw) && return false
    vals = sort(unique(all_wtp_raw), rev=true)
    length(vals) == 1 && return false

    wtp_counts = Dict{Price, Int}()
    for v in all_wtp_raw; wtp_counts[v] = get(wtp_counts, v, 0) + 1; end

    for i in 1:(length(vals) - 1)
        v      = vals[i]
        v_next = vals[i + 1]
        j_v    = wtp_counts[v]
        d_at_v = aggregate_demand(m, v)
        d_gap  = d_at_v - j_v

        s_hi = aggregate_supply(m, v)[2]         # S_max at v = min supply just above v

        # Sign change: excess demand AT v, excess supply just above v.
        # Supply just above v is S_max(v) (singleton); demand just above v is d_gap.
        # If d_gap < S_max(v): excess supply persists at ALL prices above v (supply
        # non-decreasing, demand constant in the gap).  And excess demand holds at ALL
        # prices ≤ v (demand non-decreasing downward, supply non-increasing downward).
        # So no clearing exists anywhere — bridge failure ↔ WGE fails.  No further
        # check on the gap below v is needed; the s_lo condition was a bug.
        excess_at_v    = d_at_v - s_hi
        excess_above_v = d_gap  - s_hi
        excess_at_v > 0 && excess_above_v < 0 && return true
    end
    false
end

# ── Per-size analysis ─────────────────────────────────────────────────────────

struct Obs
    j       :: Int
    bj      :: Float64
    bj2     :: Float64
    cleared :: Bool
    bf      :: Bool
end

function run_size(nc, nf; n=2_000, max_units=5, Q=100, seed_base=0)
    [begin
        rng = MersenneTwister(seed + seed_base)
        m, r = generate_good_market(rng;
            good=1, n_consumers=nc, n_firms=nf,
            max_units=max_units, Q=Q)
        Obs(jump_index(m), bj_index(m), bj2_index(m), r.cleared, bridge_failure(m))
    end for seed in 1:n]
end

function print_results(nc, nf, obs)
    n      = length(obs)
    fails  = .!([o.cleared for o in obs])
    js     = [o.j    for o in obs]
    bjs    = [o.bj   for o in obs]
    bj2s   = [o.bj2  for o in obs]
    bfs    = [o.bf   for o in obs]

    println("=" ^ 64)
    @printf("  %d agents (%dC + %dF)   n=%d\n", nc+nf, nc, nf, n)
    println("=" ^ 64)
    @printf("  P(WGE fails)           = %.4f\n", mean(fails))
    @printf("  P(bridge_failure)      = %.4f\n", mean(bfs))
    @printf("  P(BJ2 > 1)             = %.4f\n", mean(bj2s .> 1))
    @printf("  Agreement bf == fails  = %.4f\n", mean(bfs .== fails))
    @printf("  Agreement BJ2>1 == fails = %.4f\n", mean((bj2s .> 1) .== fails))
    println()

    println("  ── P(WGE fails) by predictor bin ──────────────────────────")

    # J
    println("  J:")
    for j in sort(unique(js))
        mask = js .== j
        sum(mask) < 5 && continue
        @printf("    J=%d   P(fail)=%.4f   n=%d (%.1f%%)\n",
                j, mean(fails[mask]), sum(mask), 100*mean(mask))
    end
    println()

    # BJ (relative jump vs post-jump demand), by quintile
    println("  BJ = J(v)/D(v+) [relative jump, post-jump demand]:")
    bj_qs = quantile(bjs, 0.0:0.2:1.0)
    for i in 1:5
        lo, hi = bj_qs[i], bj_qs[i+1]
        mask   = i < 5 ? (bjs .>= lo .&& bjs .< hi) : (bjs .>= lo)
        sum(mask) < 5 && continue
        @printf("    [%.3f,%.3f)  P(fail)=%.4f   n=%d\n",
                lo, hi, mean(fails[mask]), sum(mask))
    end
    println()

    # BJ2 (jump vs supply gap): threshold at 1
    println("  BJ2 = J(v)/z(v) [jump / supply gap]:")
    for (lo, hi) in [(0.0,0.5),(0.5,1.0),(1.0,2.0),(2.0,Inf)]
        mask = bj2s .>= lo .&& bj2s .< hi
        sum(mask) < 5 && continue
        hi_str = isinf(hi) ? "∞" : @sprintf("%.1f", hi)
        @printf("    [%.1f, %s)   P(fail)=%.4f   n=%d (%.1f%%)\n",
                lo, hi_str, mean(fails[mask]), sum(mask), 100*mean(mask))
    end
    println()
end

# ── Run ───────────────────────────────────────────────────────────────────────

sizes = [
    (4,   4,   40_000),
    (25,  25,  80_000),
    (50,  50, 120_000),
    (250, 250,160_000),
]

for (nc, nf, sb) in sizes
    obs = run_size(nc, nf; n=2_000, seed_base=sb)
    print_results(nc, nf, obs)
end
