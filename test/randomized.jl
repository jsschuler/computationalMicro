@testset "Randomized: find_equilibrium over 200 markets" begin
    # Non-existence is possible in multi-unit discrete markets when a demand jump
    # skips the supply level. In those cases cleared=false and price minimises |Z|.
    # The key invariant: whenever cleared=true, the price exactly clears.
    n_nonexistent = 0
    for seed in 1:200
        rng = MersenneTwister(seed)
        m, r = generate_good_market(rng;
            good=1, n_consumers=rand(rng,1:8), n_firms=rand(rng,1:8),
            max_units=rand(rng,1:6), Q=100)
        if r.cleared
            @test clears(m, r.price)
        else
            n_nonexistent += 1
            @test !clears(m, r.price)
        end
    end
    @info "Markets with no competitive equilibrium: $n_nonexistent / 200"
    @test n_nonexistent < 20
end

@testset "Randomized: tatonnement agrees with find_equilibrium" begin
    for seed in 1:200
        rng = MersenneTwister(seed)
        m, r = generate_good_market(rng;
            good=1, n_consumers=rand(rng,1:8), n_firms=rand(rng,1:8),
            max_units=rand(rng,1:6), Q=100)
        r.cleared || continue
        p_tat = tatonnement(m)
        @test clears(m, r.price)
        @test clears(m, p_tat)
    end
end

@testset "Randomized: Walras' Law over 200 single-good markets" begin
    for seed in 1:200
        rng = MersenneTwister(seed)
        m, r = generate_good_market(rng;
            good=1, n_consumers=rand(rng,1:8), n_firms=rand(rng,1:8),
            max_units=rand(rng,1:6), Q=100)
        r.cleared || continue
        mkt = Market([m], 1)
        @test walras_residual(mkt, [r.price]) == 0//1
    end
end

@testset "Randomized: multi-good Walras' Law over 100 markets" begin
    for seed in 1:100
        mkt, p_stars = generate_market(seed; k=3,
            n_consumers=rand(MersenneTwister(seed), 2:6),
            n_firms=rand(MersenneTwister(seed+1000), 2:5),
            max_units=3)
        # Walras' Law on the subset of goods with exact equilibria
        cleared = findall(r -> r.cleared, p_stars)
        isempty(cleared) && continue
        p_known = Price[p_stars[j].price for j in cleared]
        residual = sum(p_known[i] * excess_demand(mkt.goods[cleared[i]], p_known[i])
                       for i in eachindex(cleared))
        @test residual == 0//1
    end
end
