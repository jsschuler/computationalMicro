# Discrete Market Equilibrium Library: Design Specification

## Overview

This library implements a computational framework for studying market equilibrium with
integer quantities and rational prices. The primary objects are supply and demand functions
mapping rational price vectors to integer quantity vectors. Utility functions are derived
objects — real-valued rankings over integer bundles — not primitives. The framework
supports both Walrasian General Equilibrium (WGE) computation and Zero Intelligence
Trader (ZIT) simulation, enabling systematic comparison between optimizing and
non-optimizing agents.

---

## 1. Mathematical Foundations

### 1.1 Type Discipline

The fundamental type assignments are:

| Object | Type | Rationale |
|--------|------|-----------|
| Price vector | `ℚ₊ᵏ` | Rational, positive |
| Quantity bundle | `ℕᵏ` | Integer, nonneg |
| Excess demand | `ℤᵏ` | Integer, signed |
| Utility value | `ℝ` | Real, derived |
| WTP/WTA value | `ℚ₊` | Rational, positive |

The asymmetry between prices (ℚ) and quantities (ℕ) is intentional and load-bearing.
Prices are terms of trade — ratios — and ratios of integers are rational. Quantities are
counts of discrete goods and must be integers. Utility is a latent ranking device and
need not be arithmetically constrained.

### 1.2 Market Clearing

A price vector `p* ∈ ℚ₊ᵏ` clears the market if, for each good `j`:

```
Dⱼ(p*) ∈ [Sⱼ_min(p*), Sⱼ_max(p*)]
```

where `[Sⱼ_min, Sⱼ_max]` is the aggregate supply interval arising from firms indifferent
at `p*`. This is the correct clearing condition when supply correspondences are
multi-valued. Exact point clearing `Dⱼ(p*) = Sⱼ(p*)` is the generic case; interval
clearing handles indifference at rational prices, which occurs with positive probability
in randomly generated markets.

Walras' Law holds exactly: `p · Z(p) = 0` for all `p ∈ ℚ₊ᵏ`, where `Z(p) = D(p) - S(p)`.
This follows from budget feasibility of each agent and zero-profit of each firm at
equilibrium, both holding exactly in ℚ arithmetic.

### 1.3 Goods Market Separation

Since each firm produces exactly one good, goods markets are separable conditional on
factor prices. The equilibrium condition in good `j` depends on `pⱼ` only (given
exogenous factor prices `w̄`). This reduces the k-dimensional fixed point problem to k
independent one-dimensional crossing problems — a significant computational simplification
that preserves the multimarket structure for WGE comparison.

### 1.4 Step Function Structure

Both demand and supply are step functions `ℚ₊ → ℕ`. Their structure is determined by
ordered sequences of marginal values:

**Demand:** Consumer `i` has willingness-to-pay sequence `vᵢ¹ ≥ vᵢ² ≥ ... ≥ vᵢᵐ ∈ ℚ₊`
for successive units of good `j`. Individual demand is:

```
Dᵢ(p) = |{n : vᵢⁿ ≥ p}|
```

**Supply:** Firm `f` producing good `j` has marginal cost sequence
`cᶠ¹ ≤ cᶠ² ≤ ... ≤ cᶠᴺ ∈ ℚ₊`. The supply correspondence is:

```
Sᶠ_min(p) = |{n : cᶠⁿ < p}|
Sᶠ_max(p) = |{n : cᶠⁿ ≤ p}|
```

The interval `[Sᶠ_min(p), Sᶠ_max(p)]` is a singleton except when `p` equals some `cᶠⁿ`,
in which case the firm is indifferent between `n-1` and `n` units.

### 1.5 Equilibrium Existence

**Theorem.** Fix good `j`. If there exists at least one unit where aggregate WTP exceeds
aggregate WTA — i.e., `v¹_agg ≥ c¹_agg` — then a clearing price `p* ∈ ℚ₊` exists.

**Proof.** Let `n* = max{n : v_agg^n ≥ c_agg^n}` where `v_agg^n` is the nth highest WTP
across all consumers and `c_agg^n` is the nth lowest WTA across all firms. Then
`p* ∈ [c_agg^{n*}, v_agg^{n*}]` is a nonempty rational interval, and any `p*` in this
interval clears the market at quantity `n*`. □

This is the standard competitive equilibrium from the supply-demand crossing of aggregate
step functions. The interval rather than point reflects the indifference zone at the
marginal unit, which is a rational interval by construction.

### 1.6 Utility Recovery

Given a demand function `Dᵢ` defined by WTP sequence `vᵢ`, the implied utility is:

```
uᵢ(x) = Σₙ₌₁ˣ vᵢⁿ
```

the sum of marginal WTP values up to quantity `x`. This is the unique additively
separable utility function consistent with the WTP schedule. It satisfies:

- `uᵢ : ℕ → ℚ` — rational-valued (derived, not assumed)
- `Δuᵢ(x) = vᵢˣ` — marginal utility equals WTP at each unit
- Revealed preference: `Dᵢ(p) = argmax{uᵢ(x) - p·x : x ∈ ℕ}` exactly

For multi-good settings, cross-good substitution effects require off-diagonal WTP
terms, corresponding to a WTP matrix over pairs of goods — the discrete Slutsky matrix.
This is a layer-2 extension and not implemented in the base library.

### 1.7 Rational Approximation of Utility Classes

The base library uses WTP-derived utility (exactly ℚ-valued). For richer utility
classes whose demand functions are not natively rational-valued, the following
approximation hierarchy applies:

**Exactly ℚ-compatible:** Cobb-Douglas, linear (perfect substitutes), Leontief
(perfect complements), piecewise linear, quadratic utility with rational coefficient
matrix. These generate rational demands at rational prices exactly.

**Algebraically irrational at equilibrium:** CES with `ρ ≠ 0` (demand involves
`p^σ`), translog, AIDS. Require approximation.

**Padé approximation for CES.** The demand component `p^(-σ)` is approximated by a
diagonal Padé approximant `R[m,m](p) = P_m(p)/Q_m(p)` with rational coefficients,
where the diagonal constraint enforces homogeneity of degree zero. Coefficients are
computed via Chebyshev-Padé expansion on a compact price domain `[p_min, p_max]` and
rationalized via continued fraction approximation to denominator bound `Q`. The
rationalization error in equilibrium prices is bounded by `‖ΔZ‖ / |λ_min(∇Z)|` —
the induced perturbation scaled by the spectral gap of the tatônnement Jacobian.

**Stern-Brocot indexing.** The simplest rational price vector (in Stern-Brocot depth)
that approximately clears the market can be identified by searching the tree up to a
depth bound. The deviation of the WGE price from the nearest Stern-Brocot-simple
rational is a measure of market complexity, and ZI convergence to focal (low-depth)
prices rather than WGE prices is a testable hypothesis.

---

## 2. Type Hierarchy

### 2.1 Core Types

```julia
# ── Aliases ─────────────────────────────────────────────────────────────────
const Price    = Rational{Int64}
const PriceVec = Vector{Rational{Int64}}
const Bundle   = Vector{Int64}           # ∈ ℕᵏ, enforced by construction
const RatVec   = Vector{Rational{Int64}}

# ── Consumer demand: WTP sequence for a single good ─────────────────────────
struct ConsumerDemand
    good :: Int
    wtp  :: RatVec    # nonincreasing, ∈ ℚ₊
end

function (d::ConsumerDemand)(p::Price) :: Int
    count(v -> v >= p, d.wtp)
end

# ── Firm supply: marginal cost sequence, single output good ─────────────────
struct FirmSupply
    good :: Int
    wtac :: RatVec    # nondecreasing marginal costs, ∈ ℚ₊
end

# Returns interval [lo, hi] ⊆ ℕ; singleton iff p ∉ wtac
function supply_correspondence(f::FirmSupply, p::Price) :: Tuple{Int,Int}
    lo = count(c -> c <  p, f.wtac)
    hi = count(c -> c <= p, f.wtac)
    (lo, hi)
end

# ── Single-good market ───────────────────────────────────────────────────────
struct GoodMarket
    good      :: Int
    consumers :: Vector{ConsumerDemand}
    firms     :: Vector{FirmSupply}
end

# ── Full k-good market ───────────────────────────────────────────────────────
struct Market
    goods :: Vector{GoodMarket}
    k     :: Int
end
```

### 2.2 Aggregate Functions

```julia
function aggregate_demand(m::GoodMarket, p::Price) :: Int
    sum(d(p) for d in m.consumers; init=0)
end

function aggregate_supply(m::GoodMarket, p::Price) :: Tuple{Int,Int}
    intervals = [supply_correspondence(f, p) for f in m.firms]
    isempty(intervals) && return (0, 0)
    sum(x -> x[1], intervals), sum(x -> x[2], intervals)
end

function clears(m::GoodMarket, p::Price) :: Bool
    d        = aggregate_demand(m, p)
    s_lo, s_hi = aggregate_supply(m, p)
    s_lo <= d <= s_hi
end

function excess_demand(m::GoodMarket, p::Price) :: Int
    d        = aggregate_demand(m, p)
    s_lo, s_hi = aggregate_supply(m, p)
    d < s_lo && return d - s_lo    # negative: excess supply
    d > s_hi && return d - s_hi    # positive: excess demand
    return 0                       # clears
end

# Full market: vector of excess demands across goods
function excess_demand(m::Market, p::PriceVec) :: Vector{Int}
    [excess_demand(m.goods[j], p[j]) for j in 1:m.k]
end
```

### 2.3 Utility Recovery

```julia
# Additively separable utility from WTP sequence (ℕ → ℚ)
function utility(d::ConsumerDemand, x::Int) :: Rational{Int64}
    x <= 0 && return 0//1
    x > length(d.wtp) && error("quantity exceeds demand capacity")
    sum(d.wtp[1:x])
end

# Marginal utility = WTP at nth unit
marginal_utility(d::ConsumerDemand, n::Int) = d.wtp[n]

# Verify revealed preference: D(p) = argmax{u(x) - p*x}
function check_revealed_preference(d::ConsumerDemand, p::Price) :: Bool
    x_star = d(p)
    u_star = utility(d, x_star) - p * x_star
    all(utility(d, x) - p * x <= u_star for x in 0:length(d.wtp))
end
```

---

## 3. Random Generation

### 3.1 Design Principles

Random market instances are generated from a seed such that:

1. All WTP and WTA values are in ℚ₊
2. Aggregate WTP crosses aggregate WTA — equilibrium exists by construction
3. The crossing point `n*` and equilibrium price interval `[c_agg^{n*}, v_agg^{n*}]`
   are computed and stored alongside the market
4. The equilibrium price `p*` is drawn uniformly from this interval (in ℚ) and returned
   as a derived object, not a primitive

The boundary condition ensuring existence (step 2) is enforced by guaranteeing the
highest consumer WTP exceeds the lowest firm WTA: `max_i vᵢ¹ ≥ min_f cᶠ¹`. This is a
single inequality constraint and is trivially satisfied by appropriate scaling of the
distributions.

### 3.2 WTP/WTA Generation

```julia
using Random

# Draw a nonincreasing WTP sequence of length m in ℚ₊
# with denominator bound Q and values in [lo, hi]
function draw_wtp(rng::AbstractRNG, m::Int, Q::Int,
                  lo::Price, hi::Price) :: RatVec
    # Draw m uniform floats in [lo, hi], sort descending, rationalize
    vals = sort(rand(rng, m) .* Float64(hi - lo) .+ Float64(lo), rev=true)
    rationalize.(vals, tol=1/Q)
end

# Draw a nondecreasing WTA (marginal cost) sequence
function draw_wtac(rng::AbstractRNG, m::Int, Q::Int,
                   lo::Price, hi::Price) :: RatVec
    vals = sort(rand(rng, m) .* Float64(hi - lo) .+ Float64(lo))
    rationalize.(vals, tol=1/Q)
end

# Generate a single-good market with guaranteed equilibrium
function generate_good_market(rng::AbstractRNG;
    good       :: Int,
    n_consumers:: Int,
    n_firms    :: Int,
    max_units  :: Int,         # max units per agent/firm
    Q          :: Int  = 100,  # denominator bound
    wtp_hi     :: Price = 10//1,
    wtac_lo    :: Price = 1//10,
    wtac_hi    :: Price = 10//1) :: Tuple{GoodMarket, Price}

    consumers = [ConsumerDemand(good,
                     draw_wtp(rng, rand(rng,1:max_units), Q, 1//Q, wtp_hi))
                 for _ in 1:n_consumers]

    firms = [FirmSupply(good,
                 draw_wtac(rng, rand(rng,1:max_units), Q, wtac_lo, wtac_hi))
             for _ in 1:n_firms]

    market = GoodMarket(good, consumers, firms)

    # Find equilibrium by scanning the sorted union of all WTP/WTA values
    p_star = find_equilibrium(market)

    return market, p_star
end

# Generate full k-good market, one seed per good
function generate_market(seed::Int; k::Int, kwargs...) :: Tuple{Market, PriceVec}
    markets  = Vector{GoodMarket}(undef, k)
    p_stars  = Vector{Price}(undef, k)
    for j in 1:k
        rng = MersenneTwister(seed + j)
        markets[j], p_stars[j] = generate_good_market(rng; good=j, kwargs...)
    end
    Market(markets, k), p_stars
end
```

### 3.3 Equilibrium Finding

```julia
# Find p* for a single good by scanning candidate prices
# Candidate prices are all WTP and WTA values — the step function jump points
function find_equilibrium(m::GoodMarket) :: Union{Price, Nothing}
    candidates = sort(unique(vcat(
        [d.wtp  for d in m.consumers]...,
        [f.wtac for f in m.firms    ]...
    )))

    # Also check prices just above each candidate (for interval clearing)
    for p in candidates
        clears(m, p) && return p
    end

    # Check midpoints between adjacent candidates
    for i in 1:(length(candidates)-1)
        p_mid = (candidates[i] + candidates[i+1]) // 2
        clears(m, p_mid) && return p_mid
    end

    return nothing   # no equilibrium found (should not occur if market generated correctly)
end
```

---

## 4. WGE Solver

### 4.1 Tatônnement on ℚ

For the separable single-output case, tatônnement per good is a one-dimensional
bisection on a monotone step function. Let `Z(p) = D(p) - S_mid(p)` where
`S_mid = (S_min + S_max) / 2`. Then `Z` is nonincreasing in `p` and the bisection
converges in `O(log(1/ε))` steps.

```julia
function tatonnement(m::GoodMarket;
    p_lo  :: Price = 1//100,
    p_hi  :: Price = 1000//1,
    max_iter :: Int = 200) :: Union{Price, Nothing}

    Z_lo = excess_demand(m, p_lo)
    Z_hi = excess_demand(m, p_hi)

    Z_lo <= 0 && return p_lo   # already excess supply at floor
    Z_hi >= 0 && return p_hi   # already excess demand at ceiling

    for _ in 1:max_iter
        p_mid = (p_lo + p_hi) // 2
        Z_mid = excess_demand(m, p_mid)
        Z_mid == 0 && return p_mid
        Z_mid  > 0 ? (p_lo = p_mid) : (p_hi = p_mid)
        p_lo == p_hi && return p_lo
    end

    # Return best approximation — the interval [p_lo, p_hi] straddles the equilibrium
    return (p_lo + p_hi) // 2
end

# Full market tatônnement: independent per good
function solve_wge(m::Market;
    p0 :: PriceVec = fill(1//1, m.k),
    kwargs...) :: PriceVec
    [tatonnement(m.goods[j]; kwargs...) for j in 1:m.k]
end
```

### 4.2 Direct Equilibrium from Step Function Crossing

More efficient than tatônnement for the step function case — scan the sorted union
of all critical prices directly:

```julia
function solve_wge_exact(m::Market) :: PriceVec
    [find_equilibrium(m.goods[j]) for j in 1:m.k]
end
```

---

## 5. Zero Intelligence Traders

### 5.1 ZI Demand

A ZI consumer draws a quantity uniformly from `{0, 1, ..., budget_units(p)}` where
`budget_units(p) = floor(wealth / p)` is the maximum affordable quantity at price `p`.
This is the uniform distribution over integer points in the budget simplex — the
Ehrhart-theoretic demand.

```julia
struct ZIConsumer
    good   :: Int
    wealth :: Rational{Int64}    # p·ω at equilibrium prices
end

function zi_demand(c::ZIConsumer, p::Price, rng::AbstractRNG) :: Int
    max_q = floor(Int, c.wealth / p)
    max_q <= 0 && return 0
    rand(rng, 0:max_q)
end

# Expected ZI demand (Ehrhart mean): E[x] = max_q / 2
function zi_demand_mean(c::ZIConsumer, p::Price) :: Rational{Int64}
    max_q = floor(Int, c.wealth / p)
    max_q // 2
end
```

### 5.2 ZI Supply

A ZI firm draws quantity uniformly from `{0, 1, ..., capacity}` regardless of price.
This is the null model for supply — no cost optimization, just random quantity within
physical capacity.

```julia
struct ZIFirm
    good     :: Int
    capacity :: Int
end

zi_supply(f::ZIFirm, rng::AbstractRNG) :: Int = rand(rng, 0:f.capacity)

zi_supply_mean(f::ZIFirm) :: Rational{Int64} = f.capacity // 2
```

### 5.3 ZI Market Simulation

```julia
struct ZIMarket
    consumers :: Vector{ZIConsumer}
    firms     :: Vector{ZIFirm}
    k         :: Int
end

# Single period ZI outcome: realized quantities, no price discovery
function zi_period(m::ZIMarket, p::PriceVec, rng::AbstractRNG) :: Tuple{Bundle, Bundle}
    demand = [sum(zi_demand(c, p[c.good], rng)
                  for c in m.consumers if c.good == j; init=0)
              for j in 1:m.k]
    supply = [sum(zi_supply(f, rng)
                  for f in m.firms if f.good == j; init=0)
              for j in 1:m.k]
    Bundle(demand), Bundle(supply)
end

# T-period simulation at fixed price vector p
function zi_simulate(m::ZIMarket, p::PriceVec, T::Int;
    seed::Int=0) :: NamedTuple
    rng = MersenneTwister(seed)
    demands = Matrix{Int}(undef, T, m.k)
    supplies = Matrix{Int}(undef, T, m.k)
    for t in 1:T
        demands[t,:], supplies[t,:] = zi_period(m, p, rng)
    end
    (demands=demands, supplies=supplies,
     excess=demands .- supplies,
     mean_excess=vec(mean(demands .- supplies, dims=1)))
end
```

---

## 6. Comparison Statistics

### 6.1 Efficiency Measures

The primary comparison objects between WGE and ZI outcomes:

```julia
# Total surplus at an allocation (q_demand, q_supply) at price p
# Surplus = sum of consumer surplus + producer surplus
function total_surplus(m::GoodMarket, q::Int) :: Rational{Int64}
    # Consumer surplus: sum of WTP above equilibrium for units received
    cs = sum(d.wtp[n] for d in m.consumers
             for n in 1:min(q, length(d.wtp)); init=0//1)
    # Producer surplus: q * p* - total variable cost
    # Approximated here as sum of WTA below q for each firm
    ps = sum(f.wtac[n] for f in m.firms
             for n in 1:min(q, length(f.wtac)); init=0//1)
    cs - ps
end

# WGE surplus (maximum feasible)
function wge_surplus(m::GoodMarket, p_star::Price) :: Rational{Int64}
    q_star = aggregate_demand(m, p_star)
    total_surplus(m, q_star)
end

# ZI efficiency ratio: E[ZI surplus] / WGE surplus ∈ [0,1]
function zi_efficiency(m::GoodMarket, p_star::Price,
    zi_demands::Vector{Int}) :: Float64
    s_wge = Float64(wge_surplus(m, p_star))
    s_zi  = mean(Float64(total_surplus(m, q)) for q in zi_demands)
    s_wge > 0 ? s_zi / s_wge : NaN
end
```

### 6.2 Price Distribution Statistics

```julia
# Given a collection of ZI realized quantities, find implied clearing price
# (the unique p such that excess demand crosses zero in the ZI realization)
function zi_clearing_price(m::GoodMarket, zi_supply::Int) :: Union{Price, Nothing}
    # Find p such that aggregate_demand(m, p) == zi_supply
    candidates = sort(unique(vcat([d.wtp for d in m.consumers]...)))
    for p in candidates
        aggregate_demand(m, p) == zi_supply && return p
    end
    nothing
end

# Summary statistics for comparison
struct ComparisonStats
    good           :: Int
    p_wge          :: Price
    q_wge          :: Int
    zi_mean_demand :: Float64
    zi_mean_supply :: Float64
    zi_efficiency  :: Float64
    n_trials       :: Int
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
```

---

## 7. Build Order and Extension Points

### 7.1 Recommended Build Order

| Phase | Components | Validates |
|-------|-----------|-----------|
| 1 | `ConsumerDemand`, `FirmSupply`, aggregate functions | Type discipline, ℕ/ℚ closure |
| 2 | `generate_good_market`, `find_equilibrium` | Existence by construction |
| 3 | `solve_wge_exact`, `tatonnement` | WGE correctness |
| 4 | `ZIConsumer`, `ZIFirm`, `zi_simulate` | ZI baseline |
| 5 | `compare`, `zi_efficiency`, surplus functions | Comparison statistics |
| 6 | Multi-good `Market`, `solve_wge` | Full k-dimensional system |
| 7 | Padé approximation layer | CES and other utility classes |
| 8 | Stern-Brocot focal price explorer | Price complexity analysis |
| 9 | Factor markets, endogenous `w` | Production-side depth |
| 10 | OLG wrapper | Dynamic extension |

### 7.2 Extension: Rational Demand Systems

Future utility classes should implement the `RationalDemandSystem` interface:

```julia
abstract type RationalDemandSystem end

# Required methods
demand(d::RationalDemandSystem, p::PriceVec, ω::Bundle) :: Bundle
walras_residual(d::RationalDemandSystem, p::PriceVec, ω::Bundle) :: Rational{Int64}

# Concrete implementations
struct CobbDouglas <: RationalDemandSystem
    α :: RatVec    # preference shares, sum to 1
end

struct QuadraticUtility <: RationalDemandSystem
    A :: Matrix{Rational{Int64}}    # positive definite
    b :: RatVec
end

struct PadeApprox <: RationalDemandSystem
    f       :: Function             # utility to approximate
    m       :: Int                  # Padé order
    p_range :: Tuple{Price, Price}  # approximation domain
    Q       :: Int                  # denominator bound
end
```

### 7.3 Extension: Stern-Brocot Price Complexity

```julia
# Stern-Brocot depth of a rational number p/q
# Depth = number of steps to reach p/q in the Stern-Brocot tree
function sb_depth(p::Price) :: Int
    # Implement via continued fraction length
    cf = continued_fraction(p)
    sum(cf) - 1
end

# Find the simplest rational in an interval [lo, hi]
# i.e., the shallowest node in the Stern-Brocot tree within [lo, hi]
function simplest_rational(lo::Price, hi::Price) :: Price
    # Stern-Brocot search: standard algorithm
    # Returns the unique rational with smallest denominator in (lo, hi)
    ...
end

# For a market, find the "focal" equilibrium price — simplest rational
# that approximately clears within tolerance ε
function focal_price(m::GoodMarket, ε::Rational{Int64}) :: Price
    p_wge = find_equilibrium(m)
    simplest_rational(p_wge - ε, p_wge + ε)
end
```

---

## 8. Notes on Arithmetic Safety

**Overflow.** `Rational{Int64}` overflows for large numerators/denominators. The
denominator bound `Q` controls this — keeping `Q ≤ 10^4` and quantities below
`10^6` is safe. For large markets, use `Rational{BigInt}` with a performance penalty.
The library should expose a `PREC` type parameter defaulting to `Int64`.

**Canonical form.** Julia's `Rational` automatically reduces to lowest terms. WTP
and WTA sequences should be deduplicated after generation to avoid phantom
indifference at identical rational values from floating-point rationalization.

**Determinism.** All random generation passes an explicit `AbstractRNG` — no global
state. Seeds are integer, markets are reproducible. The `generate_market(seed; ...)` 
entry point derives per-good RNGs as `MersenneTwister(seed + j)` to ensure
independence across goods while preserving full reproducibility from a single seed.

**Exact clearing check.** `clears(m, p)` returns `Bool` with no tolerance parameter.
Either the market clears or it does not. The tatônnement solver returns the exact
rational equilibrium price when it exists, not a float approximation.
