@testset "Phase 1: demand / supply / aggregate" begin

    c = ConsumerDemand(1, [6//1, 4//1, 2//1])
    f = FirmSupply(1,   [1//1, 3//1, 5//1])
    m = GoodMarket(1, [c], [f])

    @test c(6//1) == 1
    @test c(4//1) == 2
    @test c(2//1) == 3
    @test c(7//1) == 0

    @test supply_correspondence(f, 3//1) == (1, 2)
    @test supply_correspondence(f, 2//1) == (1, 1)
    @test supply_correspondence(f, 1//1) == (0, 1)
    @test supply_correspondence(f, 6//1) == (3, 3)

    @test aggregate_demand(m, 4//1) == 2
    @test aggregate_supply(m, 3//1) == (1, 2)

    @test  clears(m, 3//1)
    @test  clears(m, 4//1)
    @test !clears(m, 5//1)

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

@testset "Phase 1: walras residual" begin

    # Single-good market, trivially satisfies p·Z = 0 at equilibrium
    c1 = ConsumerDemand(1, [9//1, 7//1, 5//1])
    f1 = FirmSupply(1,   [3//1, 5//1, 7//1])
    gm = GoodMarket(1, [c1], [f1])
    mkt = Market([gm], 1)

    # At any price where market clears, walras_residual == 0
    for p in [3//1, 5//1, 7//1]
        if clears(gm, p)
            @test walras_residual(mkt, [p]) == 0//1
        end
    end

end
