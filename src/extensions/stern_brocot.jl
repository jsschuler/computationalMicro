# Phase 8 (optional extension): Stern-Brocot focal price explorer
#
# Load with:
#   include(joinpath(pkgdir(DiscreteMarket), "src", "extensions", "stern_brocot.jl"))
#
# Depends on: core.jl, solvers.jl

import DiscreteMarket: find_equilibrium

# ── Continued fraction expansion ──────────────────────────────────────────────

# Standard (greedy) continued fraction of p > 0, via Euclidean algorithm.
# Result uses the canonical form where the last coefficient ≥ 1 (and ≥ 2 for p ≠ 1).
function continued_fraction(p::Price) :: Vector{Int}
    p > 0//1 || error("price must be positive")
    n, d = numerator(p), denominator(p)
    result = Int[]
    while d != 0
        a, r = divrem(n, d)
        push!(result, a)
        n, d = d, r
    end
    result
end

# Stern-Brocot tree depth: number of left/right steps to reach p.
# Equals sum(continued_fraction(p)) - 1.
# Examples: 1//1→0, 1//2→1, 2//3→2, 3//5→3.
function sb_depth(p::Price) :: Int
    sum(continued_fraction(p)) - 1
end

# ── Simplest rational in a closed interval ────────────────────────────────────

# Returns the rational in [lo, hi] with the minimum denominator,
# i.e., the shallowest Stern-Brocot node in the interval.
# Algorithm: if the interval straddles an integer, return that integer;
# otherwise shift to (0,1) via x ↦ 1/(x−n) and recurse.
function simplest_rational(lo::Price, hi::Price) :: Price
    lo <= hi || error("empty interval: lo=$lo > hi=$hi")
    lo == hi && return lo

    # Check for an integer in [lo, hi]
    a = ceil(Int, lo)
    a <= floor(Int, hi) && return Price(a)

    # lo and hi are in the same open unit interval (n, n+1)
    n = floor(Int, lo)
    # Map via x ↦ 1/(x − n): the simplest r ∈ [lo, hi] is n + 1/s
    # where s is the simplest rational in [1/(hi−n), 1/(lo−n)]
    inner = simplest_rational(inv(hi - n), inv(lo - n))
    (n * numerator(inner) + denominator(inner)) // numerator(inner)
end

# ── Focal price and complexity ────────────────────────────────────────────────

# The "focal" equilibrium price: shallowest Stern-Brocot rational that lies
# within ε of the WGE price.  Returns nothing when no equilibrium exists.
function focal_price(m, ε::Price = 1//100) :: Union{Price, Nothing}
    r = find_equilibrium(m)
    r.cleared || return nothing
    simplest_rational(r.price - ε, r.price + ε)
end

# Stern-Brocot depth of the WGE price — a measure of market "arithmetic complexity".
function price_complexity(m) :: Union{Int, Nothing}
    r = find_equilibrium(m)
    r.cleared || return nothing
    sb_depth(r.price)
end

struct FocalPriceStats
    p_wge        :: Union{Price, Nothing}
    p_focal      :: Union{Price, Nothing}
    depth_wge    :: Union{Int,   Nothing}
    depth_focal  :: Union{Int,   Nothing}
    depth_saving :: Union{Int,   Nothing}   # depth_wge − depth_focal ≥ 0
end

function focal_stats(m, ε::Price = 1//100) :: FocalPriceStats
    r = find_equilibrium(m)
    r.cleared || return FocalPriceStats(nothing, nothing, nothing, nothing, nothing)
    p_focal = simplest_rational(r.price - ε, r.price + ε)
    d_wge   = sb_depth(r.price)
    d_focal = sb_depth(p_focal)
    FocalPriceStats(r.price, p_focal, d_wge, d_focal, d_wge - d_focal)
end
