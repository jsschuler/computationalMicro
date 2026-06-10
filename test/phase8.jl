using Statistics
include(joinpath(pkgdir(DiscreteMarket), "src", "extensions", "stern_brocot.jl"))

# ── continued_fraction ────────────────────────────────────────────────────────

@testset "Phase 8: continued_fraction" begin
    @test continued_fraction(1//1) == [1]
    @test continued_fraction(2//1) == [2]
    @test continued_fraction(1//2) == [0, 2]
    @test continued_fraction(3//2) == [1, 2]
    @test continued_fraction(2//3) == [0, 1, 2]
    @test continued_fraction(5//3) == [1, 1, 2]
    @test continued_fraction(5//7) == [0, 1, 2, 2]
    @test continued_fraction(7//5) == [1, 2, 2]
    # Reconstruction: evaluate CF back to rational
    function eval_cf(cf)
        isempty(cf) && return 0//1
        r = cf[end] // 1
        for a in cf[end-1:-1:1]
            r = a + 1//r
        end
        r
    end
    for p in Price[1//1, 3//2, 7//5, 5//7, 13//8, 21//13]
        @test eval_cf(continued_fraction(p)) == p
    end
end

# ── sb_depth ──────────────────────────────────────────────────────────────────

@testset "Phase 8: sb_depth" begin
    # Depth 0: root
    @test sb_depth(1//1) == 0

    # Depth 1
    @test sb_depth(1//2) == 1
    @test sb_depth(2//1) == 1

    # Depth 2
    @test sb_depth(1//3) == 2
    @test sb_depth(2//3) == 2
    @test sb_depth(3//2) == 2
    @test sb_depth(3//1) == 2

    # Depth 3: cf sum = 4
    @test sb_depth(3//5) == 3   # cf = [0,1,1,2]
    @test sb_depth(5//3) == 3   # cf = [1,1,2]
    @test sb_depth(5//7) == 4   # cf=[0,1,2,2], sum=5, depth=4

    # Verify depth monotonicity: p//1 has depth p-1
    for n in 1:6
        @test sb_depth(n//1) == n - 1
    end

    # Verify symmetry: sb_depth(p//q) == sb_depth(q//p)
    for p in Price[2//3, 3//5, 5//7, 7//11]
        @test sb_depth(p) == sb_depth(inv(p))
    end
end

# ── simplest_rational ─────────────────────────────────────────────────────────

@testset "Phase 8: simplest_rational — basic cases" begin
    # Integer in interval → returns smallest integer
    @test simplest_rational(1//2, 3//2) == 1//1
    @test simplest_rational(3//4, 5//2) == 1//1
    @test simplest_rational(5//2, 7//2) == 3//1

    # Endpoint is the simplest
    @test simplest_rational(1//3, 1//2) == 1//2    # depth 1 < depth 2
    @test simplest_rational(2//3, 3//4) == 2//3    # depth 2 < depth 3
    @test simplest_rational(3//5, 2//3) == 2//3    # 2//3 depth 2 < 3//5 depth 3

    # Point interval
    @test simplest_rational(3//7, 3//7) == 3//7

    # Wide interval contains a simple rational
    @test simplest_rational(1//4, 3//4) == 1//2
    @test simplest_rational(2//5, 3//4) == 1//2
end

@testset "Phase 8: simplest_rational — depth is minimal" begin
    # For each test interval, verify the returned rational has depth ≤ any
    # other rational we can find in the interval.
    function check_minimal_depth(lo, hi)
        s = simplest_rational(lo, hi)
        d_s = sb_depth(s)
        # Check a sample of candidate rationals in [lo, hi]
        candidates = [n//d for d in 1:20 for n in floor(Int,lo*d):ceil(Int,hi*d)
                       if lo <= n//d <= hi && gcd(n, d) == 1 && n > 0]
        all(sb_depth(c) >= d_s for c in candidates)
    end

    @test check_minimal_depth(1//3, 1//2)
    @test check_minimal_depth(2//5, 3//4)
    @test check_minimal_depth(3//7, 5//11)
    @test check_minimal_depth(7//10, 9//10)
    @test check_minimal_depth(4//5,  6//7)
end

@testset "Phase 8: simplest_rational — focal prices are not deeper than WGE" begin
    # simplest_rational(p - ε, p + ε) should have depth ≤ sb_depth(p)
    for p in Price[3//2, 5//3, 7//5, 5//7, 8//5]
        s = simplest_rational(p - 1//20, p + 1//20)
        @test sb_depth(s) <= sb_depth(p)
    end
end

# ── focal_price and focal_stats ───────────────────────────────────────────────

@testset "Phase 8: focal_price — recovers WGE for tight ε" begin
    # p* = 3//1 for the simple market (depth 2); ε = 1//4 keeps us at 3//1
    m = GoodMarket(1,
        [ConsumerDemand(1, [6//1, 4//1, 2//1])],
        [FirmSupply(1,   [1//1, 3//1, 5//1])])
    @test focal_price(m, 1//4) == 3//1
    @test focal_price(m, 1//100) == 3//1
end

@testset "Phase 8: focal_price — finds simpler price for wide ε" begin
    # Build a market whose WGE price is 7//5 (depth 3).
    # WTP:  7//5 (one consumer, one unit)
    # WTA:  7//5 (one firm, one unit)
    # Equilibrium at 7//5.  With ε = 1//5, the interval [6//5, 8//5] contains
    # 1//1 (depth 0) and other shallower rationals.
    m = GoodMarket(1,
        [ConsumerDemand(1, [7//5])],
        [FirmSupply(1,    [7//5])])
    p_eq = find_equilibrium(m)
    @test p_eq == 7//5
    @test sb_depth(7//5) == 4   # cf=[1,2,2], depth=4

    # [6/5, 8/5] = [1.2, 1.6] — no integer; 3//2 (depth 2) lies inside
    p_f = focal_price(m, 1//5)
    @test sb_depth(p_f) < sb_depth(7//5)   # focal is shallower
    @test p_f == 3//2                        # 3//2 is the simplest in [6/5, 8/5]
end

@testset "Phase 8: focal_price — nothing when no equilibrium" begin
    # Two consumers each with WTP=4, firm with WTA=[3,5].
    # At p≤4: D=2, S_max=1 → excess demand.
    # At p∈(4,5): D=0, S_min=1 → excess supply.
    # Demand jumps from 2 to 0, skipping S=1; no equilibrium exists.
    m = GoodMarket(1,
        [ConsumerDemand(1, [4//1]), ConsumerDemand(1, [4//1])],
        [FirmSupply(1,    [3//1, 5//1])])
    @test isnothing(find_equilibrium(m))
    @test isnothing(focal_price(m))
end

@testset "Phase 8: focal_stats — depth saving is non-negative" begin
    rng = MersenneTwister(99)
    n_tested = 0
    for seed in 1:100
        rng2 = MersenneTwister(seed)
        m, p = generate_good_market(rng2;
            good=1, n_consumers=rand(rng2, 2:5), n_firms=rand(rng2, 2:5),
            max_units=3, Q=100)
        isnothing(p) && continue
        stats = focal_stats(m, 1//20)
        @test stats.depth_saving >= 0
        n_tested += 1
    end
    @test n_tested >= 50   # at least half had valid equilibria
end

@testset "Phase 8: price_complexity — sb_depth of random WGE prices" begin
    depths = Int[]
    for seed in 1:200
        rng = MersenneTwister(seed)
        m, p = generate_good_market(rng;
            good=1, n_consumers=rand(rng, 2:6), n_firms=rand(rng, 2:6),
            max_units=4, Q=100)
        isnothing(p) && continue
        push!(depths, sb_depth(p))
    end
    @test !isempty(depths)
    @info "WGE price sb_depth — min: $(minimum(depths))  mean: $(round(mean(depths), digits=2))  max: $(maximum(depths))"
    # WGE prices typically have small depth (most are simple rationals)
    @test mean(depths) < 20
end
