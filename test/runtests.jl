using Test
using Random
using DiscreteMarket

# ── helpers ──────────────────────────────────────────────────────────────────
function simple_market()
    # Consumer: WTP = [6, 4, 2]
    # Firm:     WTA = [1, 3, 5]
    # Crossing: unit 2 (WTP[2]=4 >= WTA[2]=3), unit 3 (WTP[3]=2 < WTA[3]=5)
    # => equilibrium quantity = 2, p* ∈ [3, 4]
    c = ConsumerDemand(1, [6//1, 4//1, 2//1])
    f = FirmSupply(1,   [1//1, 3//1, 5//1])
    GoodMarket(1, [c], [f])
end

# ─────────────────────────────────────────────────────────────────────────────
@testset "Phase 1: demand / supply / aggregate" begin

    m = simple_market()

    @test m.consumers[1](6//1) == 1   # only unit 1 has WTP >= 6
    @test m.consumers[1](4//1) == 2   # units 1 and 2
    @test m.consumers[1](2//1) == 3   # all three
    @test m.consumers[1](7//1) == 0   # none

    @test supply_correspondence(m.firms[1], 3//1) == (1, 2)   # c=3 is on boundary
    @test supply_correspondence(m.firms[1], 2//1) == (1, 1)   # c<2 → 1 unit strictly
    @test supply_correspondence(m.firms[1], 1//1) == (0, 1)   # c=1 on boundary
    @test supply_correspondence(m.firms[1], 6//1) == (3, 3)   # all units

    @test aggregate_demand(m, 4//1) == 2
    @test aggregate_supply(m, 3//1) == (1, 2)

    @test clears(m, 3//1)
    @test clears(m, 4//1)
    @test !clears(m, 5//1)   # D=1, S=[2,3] → excess supply

    @test excess_demand(m, 5//1) < 0
    @test excess_demand(m, 3//1) == 0
    @test excess_demand(m, 1//1) > 0

end

@testset "Phase 1: utility recovery" begin

    d = ConsumerDemand(1, [6//1, 4//1, 2//1])

    @test utility(d, 0) == 0//1
    @test utility(d, 1) == 6//1
    @test utility(d, 2) == 10//1
    @test utility(d, 3) == 12//1

    @test marginal_utility(d, 1) == 6//1
    @test marginal_utility(d, 3) == 2//1

    @test check_revealed_preference(d, 3//1)
    @test check_revealed_preference(d, 5//1)
    @test check_revealed_preference(d, 1//1)

end

@testset "Phase 2: random generation" begin

    rng = MersenneTwister(42)
    wtp = draw_wtp(rng, 5, 100, 1//10, 10//1)
    @test length(wtp) == 5
    @test all(wtp[i] >= wtp[i+1] for i in 1:4)   # nonincreasing
    @test all(v > 0 for v in wtp)

    rng2 = MersenneTwister(42)
    wtac = draw_wtac(rng2, 5, 100, 1//10, 10//1)
    @test length(wtac) == 5
    @test all(wtac[i] <= wtac[i+1] for i in 1:4)  # nondecreasing

    rng3 = MersenneTwister(7)
    gm, p = generate_good_market(rng3;
        good=1, n_consumers=3, n_firms=2, max_units=4)
    @test gm isa GoodMarket
    if !isnothing(p)
        @test clears(gm, p)
    end

end

@testset "Phase 3: equilibrium solvers" begin

    m = simple_market()

    p_exact = find_equilibrium(m)
    @test !isnothing(p_exact)
    @test clears(m, p_exact)

    p_tat = tatonnement(m)
    @test !isnothing(p_tat)
    @test clears(m, p_tat)

    mkt = Market([m], 1)
    ps_exact = solve_wge_exact(mkt)
    @test length(ps_exact) == 1
    @test !isnothing(ps_exact[1])
    @test clears(m, ps_exact[1])

    ps_tat = solve_wge(mkt)
    @test !isnothing(ps_tat[1])
    @test clears(m, ps_tat[1])

end

@testset "Phase 3: Walras' Law" begin
    # p · Z(p) == 0 at equilibrium
    m = simple_market()
    p = find_equilibrium(m)
    mkt = Market([m], 1)
    Z = excess_demand(mkt, [p])
    @test p * Z[1] == 0
end

@testset "Phase 4: ZI traders" begin

    rng = MersenneTwister(1)
    c = ZIConsumer(1, 12//1)
    f = ZIFirm(1, 6)

    q_c = zi_demand(c, 3//1, rng)      # max_q = floor(12/3) = 4
    @test 0 <= q_c <= 4

    q_f = zi_supply(f, rng)
    @test 0 <= q_f <= 6

    @test zi_demand_mean(c, 3//1) == 2//1    # 4//2
    @test zi_supply_mean(f)       == 3//1    # 6//2

    zm = ZIMarket([c], [f], 1)
    result = zi_simulate(zm, [3//1], 200; seed=0)
    @test size(result.demands)  == (200, 1)
    @test size(result.supplies) == (200, 1)
    @test size(result.excess)   == (200, 1)
    @test length(result.mean_excess) == 1

end

@testset "Phase 5: comparison statistics" begin

    m = simple_market()
    p_star = find_equilibrium(m)

    s = total_surplus(m, 2)   # q=2: CS = 6+4 = 10, PS = 1+3 = 4, surplus = 6
    @test s == 6//1

    @test wge_surplus(m, p_star) == 6//1

    zi_d = fill(2, 100)
    eff = zi_efficiency(m, p_star, zi_d)
    @test eff ≈ 1.0

    zi_d2 = fill(0, 100)
    eff2 = zi_efficiency(m, p_star, zi_d2)
    @test eff2 == 0.0

    stats = compare(m, p_star, zi_d, fill(2, 100))
    @test stats.good == 1
    @test stats.p_wge == p_star
    @test stats.q_wge == 2
    @test stats.zi_efficiency ≈ 1.0
    @test stats.n_trials == 100

end

@testset "Phase 5: multi-good market round-trip" begin

    mkt, p_stars = generate_market(42; k=3,
        n_consumers=4, n_firms=3, max_units=5)

    @test mkt.k == 3
    @test length(p_stars) == 3

    Z = excess_demand(mkt, p_stars)
    for j in 1:3
        @test clears(mkt.goods[j], p_stars[j])
        @test Z[j] == 0
    end

end

@testset "Randomized clearing: find_equilibrium over 200 markets" begin
    # Non-existence is possible in multi-unit discrete markets when a demand jump
    # skips the supply level (two consumers share a WTP, causing D to drop by 2
    # while S = 1 throughout the gap). find_equilibrium correctly returns nothing
    # in those cases. The key invariant: whenever a price IS returned, it clears.
    n_nonexistent = 0
    for seed in 1:200
        rng = MersenneTwister(seed)
        m, p = generate_good_market(rng;
            good=1, n_consumers=rand(rng,1:8), n_firms=rand(rng,1:8),
            max_units=rand(rng,1:6), Q=100)
        if isnothing(p)
            n_nonexistent += 1
        else
            @test clears(m, p)
        end
    end
    @info "Markets with no competitive equilibrium: $n_nonexistent / 200"
    # Non-existence rate should be low (well under 10%)
    @test n_nonexistent < 20
end

@testset "Randomized clearing: tatonnement agrees with find_equilibrium" begin
    disagreements = 0
    for seed in 1:200
        rng = MersenneTwister(seed)
        m, p_exact = generate_good_market(rng;
            good=1, n_consumers=rand(rng,1:8), n_firms=rand(rng,1:8),
            max_units=rand(rng,1:6), Q=100)
        isnothing(p_exact) && continue
        p_tat = tatonnement(m)
        if isnothing(p_tat)
            disagreements += 1
        else
            @test clears(m, p_tat)
            # Both must clear — they need not return the same price,
            # only the same excess demand (zero)
            @test excess_demand(m, p_exact) == 0
            @test excess_demand(m, p_tat)   == 0
        end
    end
    @test disagreements == 0
end

@testset "Randomized Walras' Law over 200 single-good markets" begin
    for seed in 1:200
        rng = MersenneTwister(seed)
        m, p = generate_good_market(rng;
            good=1, n_consumers=rand(rng,1:8), n_firms=rand(rng,1:8),
            max_units=rand(rng,1:6), Q=100)
        isnothing(p) && continue
        mkt = Market([m], 1)
        Z = excess_demand(mkt, [p])
        @test p * Z[1] == 0
    end
end
