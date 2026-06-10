@testset "Randomized: find_equilibrium over 200 markets" begin
    # Non-existence is possible in multi-unit discrete markets when a demand jump
    # skips the supply level. find_equilibrium correctly returns nothing in those
    # cases. The key invariant: whenever a price IS returned, it clears.
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
    @test n_nonexistent < 20
end

@testset "Randomized: tatonnement agrees with find_equilibrium" begin
    for seed in 1:200
        rng = MersenneTwister(seed)
        m, p_exact = generate_good_market(rng;
            good=1, n_consumers=rand(rng,1:8), n_firms=rand(rng,1:8),
            max_units=rand(rng,1:6), Q=100)
        isnothing(p_exact) && continue
        p_tat = tatonnement(m)
        @test clears(m, p_exact)
        @test clears(m, p_tat)
    end
end

@testset "Randomized: Walras' Law over 200 single-good markets" begin
    for seed in 1:200
        rng = MersenneTwister(seed)
        m, p = generate_good_market(rng;
            good=1, n_consumers=rand(rng,1:8), n_firms=rand(rng,1:8),
            max_units=rand(rng,1:6), Q=100)
        isnothing(p) && continue
        mkt = Market([m], 1)
        @test walras_residual(mkt, [p]) == 0//1
    end
end

@testset "Randomized: multi-good Walras' Law over 100 markets" begin
    for seed in 1:100
        mkt, p_stars = generate_market(seed; k=3,
            n_consumers=rand(MersenneTwister(seed), 2:6),
            n_firms=rand(MersenneTwister(seed+1000), 2:5),
            max_units=3)
        # Walras' Law on the subset of goods with known equilibria
        known = findall(!isnothing, p_stars)
        isempty(known) && continue
        p_known = Price[p_stars[j] for j in known]
        residual = sum(p_known[i] * excess_demand(mkt.goods[known[i]], p_known[i])
                       for i in eachindex(known))
        @test residual == 0//1
    end
end
