# Phase 2: random market generation
# Depends on: core.jl, solvers.jl (find_equilibrium)

function draw_wtp(rng::AbstractRNG, m::Int, Q::Int,
                  lo::Price, hi::Price) :: RatVec
    vals = sort(rand(rng, m) .* Float64(hi - lo) .+ Float64(lo), rev=true)
    rationalize.(vals, tol=1/Q)
end

function draw_wtac(rng::AbstractRNG, m::Int, Q::Int,
                   lo::Price, hi::Price) :: RatVec
    vals = sort(rand(rng, m) .* Float64(hi - lo) .+ Float64(lo))
    rationalize.(vals, tol=1/Q)
end

function generate_good_market(rng::AbstractRNG;
    good        :: Int,
    n_consumers :: Int,
    n_firms     :: Int,
    max_units   :: Int,
    Q           :: Int   = 100,
    wtp_hi      :: Price = 10//1,
    wtac_lo     :: Price = 1//10,
    wtac_hi     :: Price = 10//1) :: Tuple{GoodMarket, Union{Price,Nothing}}

    consumers = [ConsumerDemand(good,
                     draw_wtp(rng, rand(rng, 1:max_units), Q, 1//Q, wtp_hi))
                 for _ in 1:n_consumers]

    firms = [FirmSupply(good,
                 draw_wtac(rng, rand(rng, 1:max_units), Q, wtac_lo, wtac_hi))
             for _ in 1:n_firms]

    market = GoodMarket(good, consumers, firms)
    p_star = find_equilibrium(market)
    return market, p_star
end

function generate_market(seed::Int; k::Int, kwargs...) :: Tuple{Market, Vector{Union{Price,Nothing}}}
    markets = Vector{GoodMarket}(undef, k)
    p_stars = Vector{Union{Price,Nothing}}(undef, k)
    for j in 1:k
        rng = MersenneTwister(seed + j)
        m, p = generate_good_market(rng; good=j, kwargs...)
        markets[j] = m
        p_stars[j] = p
    end
    Market(markets, k), p_stars
end
