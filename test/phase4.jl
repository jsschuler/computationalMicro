@testset "Phase 4: ZI traders" begin

    rng = MersenneTwister(1)
    c = ZIConsumer(1, 12//1)
    f = ZIFirm(1, 6)

    q_c = zi_demand(c, 3//1, rng)
    @test 0 <= q_c <= 4

    q_f = zi_supply(f, rng)
    @test 0 <= q_f <= 6

    @test zi_demand_mean(c, 3//1) == 2//1
    @test zi_supply_mean(f)       == 3//1

    zm = ZIMarket([c], [f], 1)
    result = zi_simulate(zm, [3//1], 200; seed=0)
    @test size(result.demands)  == (200, 1)
    @test size(result.supplies) == (200, 1)
    @test size(result.traded)   == (200, 1)
    @test size(result.excess)   == (200, 1)
    @test length(result.mean_excess) == 1

    # traded = min(demand, supply) always
    @test all(result.traded .== min.(result.demands, result.supplies))

end
