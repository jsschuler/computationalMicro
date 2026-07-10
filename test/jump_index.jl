# Test the jump-index conjecture:
#   P(WGE exists) ≈ P(J(m) = 1)
# where J(m) = max single-price drop in aggregate demand
#            = max number of (consumer, unit) pairs sharing the same WTP value.
#
# If J = 1 everywhere, adjacent candidates always straddle equilibrium exactly
# (discrete IVT holds). If J > 1, demand can skip over any supply level in the
# gap, so WGE may not exist.
#
# Run with:
#   julia --project=. test/jump_index.jl

using Random
using Statistics
using Printf
using DiscreteMarket

# ── Jump index ────────────────────────────────────────────────────────────────
# Maximum multiplicity of any WTP value across all (consumer × unit) pairs.

function jump_index(m::GoodMarket) :: Int
    counts = Dict{Price, Int}()
    for d in m.consumers, v in d.wtp
        counts[v] = get(counts, v, 0) + 1
    end
    isempty(counts) ? 0 : maximum(values(counts))
end

# ── Per-size analysis ─────────────────────────────────────────────────────────

struct SizeResult
    n_agents   :: Int
    n_c        :: Int
    n_f        :: Int
    j_vals     :: Vector{Int}
    cleared    :: Vector{Bool}
end

function analyze_size(n_c, n_f; n_markets=1_000, max_units=5, Q=100, seed_base=0)
    j_vals  = Int[]
    cleared = Bool[]
    for seed in 1:n_markets
        rng = MersenneTwister(seed + seed_base)
        m, r = generate_good_market(rng;
            good=1, n_consumers=n_c, n_firms=n_f,
            max_units=max_units, Q=Q)
        push!(j_vals,  jump_index(m))
        push!(cleared, r.cleared)
    end
    SizeResult(n_c + n_f, n_c, n_f, j_vals, cleared)
end

# ── Print results for one size ────────────────────────────────────────────────

function print_size(sr::SizeResult)
    n   = length(sr.j_vals)
    j1  = sr.j_vals .== 1
    jgt = sr.j_vals .>  1

    p_j1       = mean(j1)
    p_wge      = mean(sr.cleared)
    p_wge_j1   = sum(sr.cleared .& j1)  / max(1, sum(j1))
    p_wge_jgt  = sum(sr.cleared .& jgt) / max(1, sum(jgt))

    println("  Agents: $(sr.n_agents)  ($(sr.n_c)C + $(sr.n_f)F)   n=$n markets")
    @printf("    P(J=1)           = %.3f\n", p_j1)
    @printf("    P(WGE exists)    = %.3f\n", p_wge)
    @printf("    P(WGE | J=1)     = %.3f\n", p_wge_j1)
    @printf("    P(WGE | J>1)     = %.3f\n", p_wge_jgt)
    @printf("    Conjecture gap   = %.4f   [P(WGE) − P(J=1)]\n", p_wge - p_j1)
    println()

    # J distribution
    max_j = min(maximum(sr.j_vals), 10)
    println("    J distribution:")
    for j in 1:max_j
        cnt  = count(==(j), sr.j_vals)
        pct  = 100 * cnt / n
        bar  = "█" ^ round(Int, 30 * cnt / n)
        @printf("      J=%2d │ %-30s %.1f%% (%d)\n", j, bar, pct, cnt)
    end
    over = count(>(max_j), sr.j_vals)
    over > 0 && @printf("      J>%2d │ %d markets\n", max_j, over)
    println()

    # Conditional WGE rate by J value (up to J=6)
    println("    P(WGE exists | J=j):")
    for j in 1:min(max_j, 6)
        mask = sr.j_vals .== j
        sum(mask) == 0 && continue
        rate = mean(sr.cleared[mask])
        @printf("      J=%d : %.3f  (n=%d)\n", j, rate, sum(mask))
    end
    println()
end

# ── Run across market sizes ───────────────────────────────────────────────────

println("=" ^ 62)
println("Jump-Index Conjecture Test")
println("  Conjecture: P(WGE exists) ≈ P(J(m) = 1)")
println("=" ^ 62)
println()

sizes = [
    (4,   4,   40_000),   # 8 agents
    (25,  25,  80_000),   # 50 agents
    (50,  50,  120_000),  # 100 agents
    (250, 250, 160_000),  # 500 agents
]

for (nc, nf, sb) in sizes
    sr = analyze_size(nc, nf; n_markets=1_000, seed_base=sb)
    print_size(sr)
end

println("=" ^ 62)
println("Interpretation guide:")
println("  Conjecture holds tightly if P(WGE) ≈ P(J=1) and")
println("  P(WGE|J=1) ≈ 1.0 while P(WGE|J>1) ≈ 0.")
println("  If P(WGE|J>1) > 0, supply can still catch demand jumps")
println("  (supply correspondence at a cost value bridges the gap).")
println("=" ^ 62)
