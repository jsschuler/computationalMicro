# Phase 3: equilibrium solvers
# Depends on: core.jl

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
    max_iter :: Int   = 200) :: Price

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

function solve_wge(m::Market; kwargs...) :: PriceVec
    [tatonnement(m.goods[j]; kwargs...) for j in 1:m.k]
end
