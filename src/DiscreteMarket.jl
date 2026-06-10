module DiscreteMarket

using Random
using Statistics

export Price, PriceVec, Bundle, RatVec
export ConsumerDemand, FirmSupply, GoodMarket, Market
export supply_correspondence, aggregate_demand, aggregate_supply, clears, excess_demand
export utility, marginal_utility, check_revealed_preference
export draw_wtp, draw_wtac, generate_good_market, generate_market
export find_equilibrium, tatonnement, solve_wge_exact, solve_wge
export ZIConsumer, ZIFirm, ZIMarket
export zi_demand, zi_demand_mean, zi_supply, zi_supply_mean
export zi_period, zi_simulate
export total_surplus, wge_surplus, zi_efficiency, zi_clearing_price
export ComparisonStats, compare

# ── Type aliases ─────────────────────────────────────────────────────────────
const Price    = Rational{Int64}
const PriceVec = Vector{Rational{Int64}}
const Bundle   = Vector{Int64}
const RatVec   = Vector{Rational{Int64}}

# ─────────────────────────────────────────────────────────────────────────────
# Phase 1: Core structs and aggregate functions
# ─────────────────────────────────────────────────────────────────────────────

struct ConsumerDemand
    good :: Int
    wtp  :: RatVec    # nonincreasing, ∈ ℚ₊
end

function (d::ConsumerDemand)(p::Price) :: Int
    count(v -> v >= p, d.wtp)
end

struct FirmSupply
    good :: Int
    wtac :: RatVec    # nondecreasing marginal costs, ∈ ℚ₊
end

function supply_correspondence(f::FirmSupply, p::Price) :: Tuple{Int,Int}
    lo = count(c -> c <  p, f.wtac)
    hi = count(c -> c <= p, f.wtac)
    (lo, hi)
end

struct GoodMarket
    good      :: Int
    consumers :: Vector{ConsumerDemand}
    firms     :: Vector{FirmSupply}
end

struct Market
    goods :: Vector{GoodMarket}
    k     :: Int
end

function aggregate_demand(m::GoodMarket, p::Price) :: Int
    sum(d(p) for d in m.consumers; init=0)
end

function aggregate_supply(m::GoodMarket, p::Price) :: Tuple{Int,Int}
    intervals = [supply_correspondence(f, p) for f in m.firms]
    isempty(intervals) && return (0, 0)
    sum(x -> x[1], intervals), sum(x -> x[2], intervals)
end

function clears(m::GoodMarket, p::Price) :: Bool
    d = aggregate_demand(m, p)
    s_lo, s_hi = aggregate_supply(m, p)
    s_lo <= d <= s_hi
end

function excess_demand(m::GoodMarket, p::Price) :: Int
    d = aggregate_demand(m, p)
    s_lo, s_hi = aggregate_supply(m, p)
    d < s_lo && return d - s_lo
    d > s_hi && return d - s_hi
    return 0
end

function excess_demand(m::Market, p::PriceVec) :: Vector{Int}
    [excess_demand(m.goods[j], p[j]) for j in 1:m.k]
end

# ─────────────────────────────────────────────────────────────────────────────
# Utility recovery
# ─────────────────────────────────────────────────────────────────────────────

function utility(d::ConsumerDemand, x::Int) :: Rational{Int64}
    x <= 0 && return 0//1
    x > length(d.wtp) && error("quantity exceeds demand capacity")
    sum(d.wtp[1:x])
end

marginal_utility(d::ConsumerDemand, n::Int) = d.wtp[n]

function check_revealed_preference(d::ConsumerDemand, p::Price) :: Bool
    x_star = d(p)
    u_star = utility(d, x_star) - p * x_star
    all(utility(d, x) - p * x <= u_star for x in 0:length(d.wtp))
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 2: Random generation
# ─────────────────────────────────────────────────────────────────────────────

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

function generate_market(seed::Int; k::Int, kwargs...) :: Tuple{Market, PriceVec}
    markets = Vector{GoodMarket}(undef, k)
    p_stars = Vector{Price}(undef, k)
    for j in 1:k
        rng = MersenneTwister(seed + j)
        m, p = generate_good_market(rng; good=j, kwargs...)
        markets[j] = m
        # Fall back to midpoint of all candidates if find_equilibrium returns nothing
        p_stars[j] = something(p, 1//1)
    end
    Market(markets, k), p_stars
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 3: Equilibrium solvers
# ─────────────────────────────────────────────────────────────────────────────

function find_equilibrium(m::GoodMarket) :: Union{Price, Nothing}
    all_vals = vcat(
        [d.wtp  for d in m.consumers]...,
        [f.wtac for f in m.firms    ]...
    )
    isempty(all_vals) && return nothing
    candidates = sort(unique(all_vals))

    for p in candidates
        clears(m, p) && return p
    end

    for i in 1:(length(candidates)-1)
        p_mid = (candidates[i] + candidates[i+1]) // 2
        clears(m, p_mid) && return p_mid
    end

    return nothing
end

function tatonnement(m::GoodMarket;
    p_lo     :: Price = 1//100,
    p_hi     :: Price = 1000//1,
    max_iter :: Int   = 200) :: Union{Price, Nothing}

    Z_lo = excess_demand(m, p_lo)
    Z_hi = excess_demand(m, p_hi)

    Z_lo <= 0 && return p_lo
    Z_hi >= 0 && return p_hi

    for _ in 1:max_iter
        # Stern-Brocot mediant: denominators grow O(n) not O(2^n), avoiding overflow
        p_mid = (numerator(p_lo) + numerator(p_hi)) //
                (denominator(p_lo) + denominator(p_hi))
        Z_mid = excess_demand(m, p_mid)
        Z_mid == 0 && return p_mid
        Z_mid  > 0 ? (p_lo = p_mid) : (p_hi = p_mid)
        p_lo == p_hi && return p_lo
    end

    return (numerator(p_lo) + numerator(p_hi)) //
           (denominator(p_lo) + denominator(p_hi))
end

function solve_wge_exact(m::Market) :: Vector{Union{Price,Nothing}}
    [find_equilibrium(m.goods[j]) for j in 1:m.k]
end

function solve_wge(m::Market; kwargs...) :: Vector{Union{Price,Nothing}}
    [tatonnement(m.goods[j]; kwargs...) for j in 1:m.k]
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 4: Zero Intelligence traders
# ─────────────────────────────────────────────────────────────────────────────

struct ZIConsumer
    good   :: Int
    wealth :: Rational{Int64}
end

function zi_demand(c::ZIConsumer, p::Price, rng::AbstractRNG) :: Int
    max_q = floor(Int, c.wealth / p)
    max_q <= 0 && return 0
    rand(rng, 0:max_q)
end

function zi_demand_mean(c::ZIConsumer, p::Price) :: Rational{Int64}
    max_q = floor(Int, c.wealth / p)
    max_q // 2
end

struct ZIFirm
    good     :: Int
    capacity :: Int
end

zi_supply(f::ZIFirm, rng::AbstractRNG) :: Int = rand(rng, 0:f.capacity)

zi_supply_mean(f::ZIFirm) :: Rational{Int64} = f.capacity // 2

struct ZIMarket
    consumers :: Vector{ZIConsumer}
    firms     :: Vector{ZIFirm}
    k         :: Int
end

function zi_period(m::ZIMarket, p::PriceVec, rng::AbstractRNG) :: Tuple{Bundle,Bundle}
    demand = [sum(zi_demand(c, p[c.good], rng)
                  for c in m.consumers if c.good == j; init=0)
              for j in 1:m.k]
    supply = [sum(zi_supply(f, rng)
                  for f in m.firms if f.good == j; init=0)
              for j in 1:m.k]
    Bundle(demand), Bundle(supply)
end

function zi_simulate(m::ZIMarket, p::PriceVec, T::Int; seed::Int=0) :: NamedTuple
    rng = MersenneTwister(seed)
    demands  = Matrix{Int}(undef, T, m.k)
    supplies = Matrix{Int}(undef, T, m.k)
    for t in 1:T
        demands[t,:], supplies[t,:] = zi_period(m, p, rng)
    end
    (demands    = demands,
     supplies   = supplies,
     excess     = demands .- supplies,
     mean_excess = vec(mean(demands .- supplies, dims=1)))
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 5: Comparison statistics
# ─────────────────────────────────────────────────────────────────────────────

function total_surplus(m::GoodMarket, q::Int) :: Rational{Int64}
    q <= 0 && return 0//1
    # Efficient allocation: assign units to highest-WTP consumers and lowest-WTA firms
    all_wtp  = sort(vcat([d.wtp  for d in m.consumers]...), rev=true)
    all_wtac = sort(vcat([f.wtac for f in m.firms    ]...))
    n = min(q, length(all_wtp), length(all_wtac))
    n == 0 && return 0//1
    sum(all_wtp[1:n]) - sum(all_wtac[1:n])
end

function wge_surplus(m::GoodMarket, p_star::Price) :: Rational{Int64}
    q_star = aggregate_demand(m, p_star)
    total_surplus(m, q_star)
end

function zi_efficiency(m::GoodMarket, p_star::Price,
                       zi_demands::Vector{Int}) :: Float64
    s_wge = Float64(wge_surplus(m, p_star))
    s_zi  = mean(Float64(total_surplus(m, q)) for q in zi_demands)
    s_wge > 0 ? s_zi / s_wge : NaN
end

function zi_clearing_price(m::GoodMarket, zi_supply::Int) :: Union{Price, Nothing}
    candidates = sort(unique(vcat([d.wtp for d in m.consumers]...)))
    for p in candidates
        aggregate_demand(m, p) == zi_supply && return p
    end
    nothing
end

struct ComparisonStats
    good            :: Int
    p_wge           :: Price
    q_wge           :: Int
    zi_mean_demand  :: Float64
    zi_mean_supply  :: Float64
    zi_efficiency   :: Float64
    n_trials        :: Int
end

function compare(m::GoodMarket, p_wge::Price,
                 zi_d::Vector{Int}, zi_s::Vector{Int}) :: ComparisonStats
    ComparisonStats(
        m.good, p_wge,
        aggregate_demand(m, p_wge),
        mean(zi_d), mean(zi_s),
        zi_efficiency(m, p_wge, zi_d),
        length(zi_d)
    )
end

end # module DiscreteMarket
