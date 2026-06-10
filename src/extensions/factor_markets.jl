# Phase 9 (optional extension): factor markets with endogenous wage
#
# Load with:
#   include(joinpath(pkgdir(DiscreteMarket), "src", "extensions", "factor_markets.jl"))
#
# Depends on: core.jl
#
# Each firm f producing good j uses labor ℓ to produce output.
# Consumer i has labor endowment eᵢ and earns wage income w·eᵢ.
# Goods markets clear at (p*, w*) jointly; factor market clears when
# labor demand from firms equals labor supply from consumers.

import DiscreteMarket: supply_correspondence, aggregate_demand, find_equilibrium

# ── Factor market types ───────────────────────────────────────────────────────

struct LaborMarket
    n_workers :: Int
    endowments :: Vector{Rational{Int64}}   # labor endowment per consumer, ∈ ℚ₊
end

# Firm with labor input: produces good `good`, uses `labor_per_unit` labor per unit
struct FirmWithLabor
    good           :: Int
    labor_per_unit :: Rational{Int64}   # ∈ ℚ₊, constant marginal labor requirement
    capacity       :: Int               # maximum units producible
end

# At wage w, firm f's marginal cost for unit n is w * labor_per_unit
function marginal_cost(f::FirmWithLabor, w::Price) :: Price
    w * f.labor_per_unit
end

# Firm supply correspondence given wage w
function supply_correspondence(f::FirmWithLabor, p::Price, w::Price) :: Tuple{Int,Int}
    mc = marginal_cost(f, w)
    lo = p > mc ? f.capacity : 0
    hi = p >= mc ? f.capacity : 0
    (lo, hi)
end

# Consumer income at wage w given labor endowment
function income(endowment::Rational{Int64}, w::Price) :: Rational{Int64}
    endowment * w
end

# ── Factor market equilibrium ─────────────────────────────────────────────────

# Aggregate labor demand from all firms at output price p and wage w
function labor_demand(firms::Vector{FirmWithLabor}, p::Price, w::Price) :: Rational{Int64}
    sum(f.labor_per_unit * supply_correspondence(f, p, w)[2]
        for f in firms; init=0//1)
end

# Labor supply = total endowment (inelastic for now)
function labor_supply(lm::LaborMarket) :: Rational{Int64}
    sum(lm.endowments)
end

# Excess labor demand at (p, w)
function labor_excess_demand(firms::Vector{FirmWithLabor},
                             lm::LaborMarket,
                             p::Price, w::Price) :: Rational{Int64}
    labor_demand(firms, p, w) - labor_supply(lm)
end

# ── Joint (p, w) equilibrium solver ──────────────────────────────────────────
#
# Algorithm: nested Stern-Brocot tatônnement.
#   Outer loop: bisect on wage w to clear the labor market.
#   Inner step: at each candidate w, convert each FirmWithLabor to an equivalent
#     FirmSupply (constant mc = w · ℓ) and solve the k goods markets independently.
#   Labor demand at (p*_j, w): the quantity actually traded in good j is
#     aggregate_demand(market_j_at_w, p*_j); the labor used by each labor firm
#     is ℓ_f · min(capacity_f, q*_j).
#
# Separability: since every firm's mc depends only on w (not on other prices),
# each goods market decouples conditional on w.  This reduces the joint fixed-point
# to k independent one-dimensional problems plus one outer bisection.

function solve_wge_with_labor(
        m::Market,
        firms_with_labor::Vector{Vector{FirmWithLabor}},
        lm::LaborMarket;
        w_lo     :: Price = 1//100,
        w_hi     :: Price = 100//1,
        max_iter :: Int   = 200) :: Tuple{Vector{EquilibriumResult}, Price}

    length(firms_with_labor) == m.k ||
        error("firms_with_labor must have one entry per good (got $(length(firms_with_labor)), need $(m.k))")

    ls = labor_supply(lm)

    # Build GoodMarket for good j at a given wage, appending labor-firm equivalents
    function market_at_w(j::Int, w::Price) :: GoodMarket
        labor_fs = [FirmSupply(j, fill(marginal_cost(f, w), f.capacity))
                    for f in firms_with_labor[j]]
        GoodMarket(j, m.goods[j].consumers, vcat(m.goods[j].firms, labor_fs))
    end

    # Solve all goods markets at wage w
    function equilibria_at(w::Price) :: Vector{EquilibriumResult}
        [find_equilibrium(market_at_w(j, w)) for j in 1:m.k]
    end

    # Labor actually used given goods equilibria at wage w.
    # Regular firms (in m.goods[j].firms) supply first (up to their hi); labor
    # firms supply the residual, up to their capacities.
    function calc_labor_demand(p_stars::Vector{EquilibriumResult}, w::Price) :: Rational{Int64}
        total = 0//1
        for j in 1:m.k
            p_j     = p_stars[j].price
            q_total = aggregate_demand(market_at_w(j, w), p_j)
            s_reg   = isempty(m.goods[j].firms) ? 0 :
                      sum(supply_correspondence(f, p_j)[2] for f in m.goods[j].firms)
            remaining = max(0, q_total - s_reg)
            for f in firms_with_labor[j]
                contrib    = min(f.capacity, remaining)
                total     += f.labor_per_unit * contrib
                remaining -= contrib
            end
        end
        total
    end

    # Tatônnement: Stern-Brocot mediant bisection on wage
    for _ in 1:max_iter
        w_mid   = (numerator(w_lo) + numerator(w_hi)) //
                  (denominator(w_lo) + denominator(w_hi))
        p_stars = equilibria_at(w_mid)
        ld      = calc_labor_demand(p_stars, w_mid)
        ez      = ld - ls

        ez == 0//1 && return (p_stars, w_mid)

        if ez > 0//1
            w_lo = w_mid   # labor too cheap → raise wage
        else
            w_hi = w_mid   # labor too expensive → lower wage
        end

        w_lo == w_hi && break
    end

    w_final = (numerator(w_lo) + numerator(w_hi)) //
              (denominator(w_lo) + denominator(w_hi))
    (equilibria_at(w_final), w_final)
end
