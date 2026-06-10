# Phase 4: Zero Intelligence traders
# Depends on: core.jl

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
    traded = min.(demands, supplies)
    (demands     = demands,
     supplies    = supplies,
     traded      = traded,
     excess      = demands .- supplies,
     mean_excess = vec(mean(demands .- supplies, dims=1)))
end
