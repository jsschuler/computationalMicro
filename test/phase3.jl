function simple_market()
    c = ConsumerDemand(1, [6//1, 4//1, 2//1])
    f = FirmSupply(1,   [1//1, 3//1, 5//1])
    GoodMarket(1, [c], [f])
end

@testset "Phase 3: equilibrium solvers" begin

    m = simple_market()

    p_exact = find_equilibrium(m)
    @test !isnothing(p_exact)
    @test clears(m, p_exact)

    p_tat = tatonnement(m)
    @test clears(m, p_tat)

    mkt = Market([m], 1)
    ps_exact = solve_wge_exact(mkt)
    @test length(ps_exact) == 1
    @test !isnothing(ps_exact[1])
    @test clears(m, ps_exact[1])

    ps_tat = solve_wge(mkt)
    @test clears(m, ps_tat[1])

end

@testset "Phase 3: Walras' Law (exact market)" begin

    m   = simple_market()
    p   = find_equilibrium(m)
    mkt = Market([m], 1)
    @test walras_residual(mkt, [p]) == 0//1

end
