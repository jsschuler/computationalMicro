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

# ── Placeholder: joint (p, w) equilibrium solver ─────────────────────────────
# Full implementation: Phase 9 — iterate tatônnement over goods prices
# and wage simultaneously, using separability conditional on w.

function solve_wge_with_labor(
        m::Market,
        firms_with_labor::Vector{Vector{FirmWithLabor}},
        lm::LaborMarket;
        w_lo :: Price = 1//100,
        w_hi :: Price = 100//1,
        max_iter :: Int = 200) :: Tuple{PriceVec, Price}
    error("Phase 9 not yet implemented: solve_wge_with_labor")
end
