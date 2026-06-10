@testset "Phase 2: random generation" begin

    rng = MersenneTwister(42)
    wtp = draw_wtp(rng, 5, 100, 1//10, 10//1)
    @test length(wtp) == 5
    @test all(wtp[i] >= wtp[i+1] for i in 1:4)
    @test all(v > 0 for v in wtp)

    rng2 = MersenneTwister(42)
    wtac = draw_wtac(rng2, 5, 100, 1//10, 10//1)
    @test length(wtac) == 5
    @test all(wtac[i] <= wtac[i+1] for i in 1:4)

    rng3 = MersenneTwister(7)
    gm, r = generate_good_market(rng3;
        good=1, n_consumers=3, n_firms=2, max_units=4)
    @test gm isa GoodMarket
    @test r isa EquilibriumResult
    if r.cleared
        @test clears(gm, r.price)
    else
        @test !clears(gm, r.price)   # min-|Z| price does not clear
    end

end
