# Phase 1: types, aggregate functions, utility recovery
# No external dependencies — loads standalone.

# ── Type aliases ──────────────────────────────────────────────────────────────
const Price    = Rational{Int64}
const PriceVec = Vector{Rational{Int64}}
const Bundle   = Vector{Int64}
const RatVec   = Vector{Rational{Int64}}

# ── Consumer demand ───────────────────────────────────────────────────────────
struct ConsumerDemand
    good :: Int
    wtp  :: RatVec    # nonincreasing, ∈ ℚ₊
end

function (d::ConsumerDemand)(p::Price) :: Int
    count(v -> v >= p, d.wtp)
end

# ── Firm supply ───────────────────────────────────────────────────────────────
struct FirmSupply
    good :: Int
    wtac :: RatVec    # nondecreasing marginal costs, ∈ ℚ₊
end

function supply_correspondence(f::FirmSupply, p::Price) :: Tuple{Int,Int}
    lo = count(c -> c <  p, f.wtac)
    hi = count(c -> c <= p, f.wtac)
    (lo, hi)
end

# ── Markets ───────────────────────────────────────────────────────────────────
struct GoodMarket
    good      :: Int
    consumers :: Vector{ConsumerDemand}
    firms     :: Vector{FirmSupply}
end

struct Market
    goods :: Vector{GoodMarket}
    k     :: Int
end

# ── Aggregate functions ───────────────────────────────────────────────────────
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

# p · Z(p) = 0 at equilibrium; exact in ℚ arithmetic
function walras_residual(m::Market, p::PriceVec) :: Rational{Int64}
    Z = excess_demand(m, p)
    sum(p[j] * Z[j] for j in 1:m.k)
end

# ── Utility recovery ──────────────────────────────────────────────────────────
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
