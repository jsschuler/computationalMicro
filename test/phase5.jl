@testset "Phase 5: comparison statistics" begin

    c = ConsumerDemand(1, [6//1, 4//1, 2//1])
    f = FirmSupply(1,   [1//1, 3//1, 5//1])
    m = GoodMarket(1, [c], [f])
    p_star = find_equilibrium(m)

    # q=2: top-2 WTP = [6,4], bottom-2 WTA = [1,3] → surplus = 10 - 4 = 6
    @test total_surplus(m, 2) == 6//1
    @test wge_surplus(m, p_star) == 6//1

    zi_d = fill(2, 100)
    zi_s = fill(2, 100)
    eff = zi_efficiency(m, p_star, zi_d)    # traded = min(d,s) = 2
    @test eff ≈ 1.0

    zi_d2 = fill(0, 100)
    eff2 = zi_efficiency(m, p_star, zi_d2)
    @test eff2 == 0.0

    stats = compare(m, p_star, zi_d, zi_s)
    @test stats.good == 1
    @test stats.p_wge == p_star
    @test stats.q_wge == 2
    @test stats.zi_efficiency ≈ 1.0
    @test stats.n_trials == 100

end

@testset "Phase 6: multi-good market" begin

    mkt, p_stars = generate_market(42; k=3,
        n_consumers=4, n_firms=3, max_units=5)

    @test mkt.k == 3
    @test length(p_stars) == 3

    stats = multi_compare(mkt, p_stars)
    @test stats.k == 3
    @test length(stats.q_wge) == 3
    @test length(stats.surplus_wge) == 3

    # For goods with a valid equilibrium: market clears and Walras holds
    for j in 1:3
        p = p_stars[j]
        isnothing(p) && continue
        @test clears(mkt.goods[j], p)
        @test excess_demand(mkt.goods[j], p) == 0
    end
    @test stats.walras == 0//1

end

@testset "Phase 6: multi-good ZI simulation" begin

    mkt, p_stars = generate_market(1; k=2,
        n_consumers=3, n_firms=2, max_units=4)

    # Only simulate on goods where equilibrium exists
    valid = findall(!isnothing, p_stars)
    @test !isempty(valid)

    p_sim = [something(p_stars[j], 1//1) for j in 1:2]

    consumers = [ZIConsumer(j, 20//1) for j in 1:2 for _ in 1:3]
    firms     = [ZIFirm(j, 5)         for j in 1:2 for _ in 1:2]
    zm = ZIMarket(consumers, firms, 2)

    result = zi_simulate(zm, p_sim, 500; seed=7)
    @test size(result.demands)  == (500, 2)
    @test size(result.traded)   == (500, 2)
    @test all(result.traded .<= result.demands)
    @test all(result.traded .<= result.supplies)

end
