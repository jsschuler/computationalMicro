# DiscreteMarket.jl

A computational framework for studying market equilibrium with integer quantities and rational prices. Supply and demand functions map rational price vectors to integer quantity vectors. Utility is a derived object, not a primitive.

## Mathematical Foundations

The type discipline is load-bearing:

| Object | Type | Rationale |
|--------|------|-----------|
| Price | `ℚ₊` | Prices are ratios — ratios of integers are rational |
| Quantity | `ℕ` | Goods are counts — counts are integers |
| Excess demand | `ℤ` | Signed integer difference |
| Utility | `ℝ` | Latent ranking device, derived from WTP |

Both demand and supply are step functions `ℚ₊ → ℕ` defined by marginal value sequences. A consumer with WTP sequence `v¹ ≥ v² ≥ ... ≥ vᵐ` demands `|{n : vⁿ ≥ p}|` units at price `p`. A firm with marginal cost sequence `c¹ ≤ c² ≤ ... ≤ cᴺ` has a supply *correspondence* `[|{n : cⁿ < p}|, |{n : cⁿ ≤ p}|]` — a singleton except when `p` hits a cost exactly.

Market clearing is interval-based: `D(p*) ∈ [S_min(p*), S_max(p*)]`. Exact point clearing is the generic case; interval clearing handles the indifference zone at rational prices. Walras' Law holds exactly: `p · Z(p) = 0` for all `p ∈ ℚ₊ᵏ`.

---

## Worked Example

### Market setup

Two consumers and two firms, each with up to 4–5 units:

```julia
c1 = ConsumerDemand(1, [9//1, 7//1, 5//1, 3//1, 1//1])  # WTP: 9,7,5,3,1
c2 = ConsumerDemand(1, [8//1, 6//1, 4//1, 2//1])          # WTP: 8,6,4,2
f1 = FirmSupply(1,   [2//1, 4//1, 6//1, 8//1])            # WTA: 2,4,6,8
f2 = FirmSupply(1,   [1//1, 3//1, 5//1, 7//1])            # WTA: 1,3,5,7
m  = GoodMarket(1, [c1, c2], [f1, f2])
```

### Supply and demand curves

Each curve is the inverse of the aggregate step function: the marginal WTP of the `q`th demanded unit (demand) and the marginal WTA of the `q`th supplied unit (supply).

```
marginal
 value
   9 │ D · · · · · · · ·
   8 │ · D · · · · · · ·
   7 │ · · D · · · · · ·
   6 │ · · · D · · · · ·
   5 │ · · · · ★ · · · ·  ← p* = 5, q* = 5
   4 │ · · · S · D · · ·
   3 │ · · S · · · D · ·
   2 │ · S · · · · · D ·
   1 │ S · · · · · · · D
     └─────────────────── quantity
       1 2 3 4 5 6 7 8 9

 D = marginal willingness-to-pay (consumers, sorted ↓)
 S = marginal willingness-to-accept (firms, sorted ↑)
 ★ = WGE: both equal 5 at q = 5
```

The `★` marks where the two sequences cross: the 5th-highest WTP (= 5) equals the 5th-lowest WTA (= 5). Every unit to the left of `★` is a profitable trade; every unit to the right destroys value.

### Aggregate step functions

```julia
p_star = find_equilibrium(m)   # => 5//1
```

| Price | D(p) | S_min(p) | S_max(p) | Z(p) |      |
|------:|-----:|---------:|---------:|-----:|------|
|     1 |    9 |        0 |        1 |   +8 |      |
|     2 |    8 |        1 |        2 |   +6 |      |
|     3 |    7 |        2 |        3 |   +4 |      |
|     4 |    6 |        3 |        4 |   +2 |      |
| **5** |**5** |      **4**|      **5**|  **0**| **p\*** |
|     6 |    4 |        5 |        6 |   −1 |      |
|     7 |    3 |        6 |        7 |   −3 |      |
|     8 |    2 |        7 |        8 |   −5 |      |
|     9 |    1 |        8 |        8 |   −7 |      |

`Z` is monotone nonincreasing; the crossing is exact at `p* = 5`.

### WGE solution

```
p* = 5,   q* = 5,   total surplus = 20
  consumer surplus = 10   (= Σ WTP for units bought − p*·q*)
  producer surplus = 10   (= p*·q* − Σ WTA for units sold)
```

Surplus by quantity (efficient allocation of `q` units across agents):

```
  q=0  │                    surplus = 0
  q=1  │ ████████            surplus = 8   (WTP 9 − WTA 1)
  q=2  │ ██████████████      surplus = 14
  q=3  │ ██████████████████  surplus = 18
  q=4  │ ████████████████████surplus = 20  ← max
  q=5  │ ████████████████████surplus = 20  ← p* (5th WTP = 5th WTA)
  q=6  │ ██████████████████  surplus = 18
  q=7  │ ██████████████      surplus = 14
  q=8  │ ████████            surplus = 8
       └──────────────────────────────────
         0       10       20
```

The plateau at `q = 4–5` reflects the marginal unit being exactly break-even (WTP = WTA = 5).

---

## Zero Intelligence Traders

ZI consumers draw quantity uniformly from `{0, 1, ..., floor(wealth/p)}`. ZI firms draw uniformly from `{0, 1, ..., capacity}`. Neither optimizes; both ignore prices except as budget constraints.

```julia
zm = ZIMarket(
    [ZIConsumer(1, 25//1), ZIConsumer(1, 25//1)],  # wealth = 25 each → max 5 units at p*=5
    [ZIFirm(1, 5),         ZIFirm(1, 5)],           # capacity = 5 each
    1)

result = zi_simulate(zm, [5//1], 2000; seed=42)
# result.mean_excess[1] ≈ 0.07   (unbiased around zero — ZI demand ≈ ZI supply on average)
```

### ZI outcome distribution (2000 trials at p* = 5)

Realized trades = min(ZI demand, ZI supply) per period:

```
 q= 0 │ ███                    124 trials
 q= 1 │ █████                  199
 q= 2 │ ███████                271
 q= 3 │ █████████              360
 q= 4 │ █████████              398  ← mode
 q= 5 │ ███████                306  ← WGE quantity
 q= 6 │ █████                  175
 q= 7 │ ██                     109
 q= 8 │ █                       46
 q= 9 │ ·                       12
       └────────────────────────────
```

### WGE vs ZI comparison

| Metric | WGE | ZI (mean over 2000 trials) |
|--------|----:|---------------------------:|
| Quantity | 5 | 4.1 (realized trades) |
| Surplus | 20 | 15.5 |
| Efficiency | 100% | **77.7%** |

ZI achieves about 78% of the efficient surplus. The loss comes from dispersion: roughly 16% of periods trade fewer than 4 units (stranding profitable trades) and another 8% trade more than 5 (executing loss-making units).

```julia
zi_eff = zi_efficiency(m, 5//1, min.(result.demands[:,1], result.supplies[:,1]))
# => 0.777
```

---

## Implementation Status

Phases 1–8 of the [design spec](discrete_market_spec.md) are complete and tested. Phases 7–8 are optional extensions loaded via `include`.

### Phase 1 — Core Types and Aggregate Functions

```julia
c = ConsumerDemand(1, [6//1, 4//1, 2//1])
c(4//1)                          # => 2  (units with WTP >= 4)

f = FirmSupply(1, [1//1, 3//1, 5//1])
supply_correspondence(f, 3//1)   # => (1, 2)  (indifferent at c=3)

m = GoodMarket(1, [c], [f])
aggregate_demand(m, 4//1)        # => 2
aggregate_supply(m, 3//1)        # => (1, 2)
clears(m, 3//1)                  # => true
excess_demand(m, 5//1)           # => -1
```

Utility is recovered from the WTP sequence as `u(x) = Σₙ₌₁ˣ vⁿ`. Revealed preference holds exactly:

```julia
utility(c, 2)                      # => 10//1
check_revealed_preference(c, 3//1) # => true
```

### Phase 2 — Random Market Generation

```julia
rng = MersenneTwister(42)
market, p_star = generate_good_market(rng;
    good=1, n_consumers=5, n_firms=3, max_units=4, Q=100)
# p_star :: Union{Price, Nothing} — nothing if no competitive equilibrium exists

mkt, p_stars = generate_market(42; k=3, n_consumers=4, n_firms=3, max_units=5)
# p_stars :: Vector{Union{Price, Nothing}} — propagates nothing per good
```

### Phase 3 — WGE Solvers

**Exact scan** (`find_equilibrium`) — scans the sorted union of all WTP/WTA values. Returns `Union{Price, Nothing}`: `nothing` when no competitive equilibrium exists. `O(n log n)`.

**Tatônnement** (`tatonnement`) — bisection on ℚ using the Stern-Brocot mediant to keep denominators bounded. Always returns a `Price` — the exact clearing price when found, or the best mediant approximation after `max_iter` steps. `solve_wge` returns `PriceVec`.

```julia
find_equilibrium(m)    # :: Union{Price, Nothing} — exact or nothing
tatonnement(m)         # :: Price                 — always returns a price
solve_wge_exact(mkt)   # :: Vector{Union{Price, Nothing}}
solve_wge(mkt)         # :: PriceVec
```

### Phase 4 — Zero Intelligence Traders

```julia
c_zi = ZIConsumer(1, 12//1)   # wealth = 12
f_zi = ZIFirm(1, 6)           # capacity = 6

zi_demand_mean(c_zi, 3//1)    # => 2//1  (floor(12/3) / 2)
zi_supply_mean(f_zi)          # => 3//1  (6 / 2)

zm = ZIMarket([c_zi], [f_zi], 1)
result = zi_simulate(zm, [3//1], 1000; seed=0)
```

### Phase 5 — Comparison Statistics

`total_surplus(m, q)` ranks all WTP values descending and all WTA values ascending, then computes `Σ WTP[1:n] - Σ WTA[1:n]` where `n = min(q, total_wtp_units, total_wta_units)` — the efficient allocation of `q` units.

```julia
total_surplus(m, 5)    # => 20//1  (exact rational)
wge_surplus(m, 5//1)   # => 20//1  (calls total_surplus at q*)

# zi_efficiency takes realized trades: min(demand, supply) per period
traded = min.(result.demands[:,1], result.supplies[:,1])
zi_efficiency(m, 5//1, traded)   # => 0.777

# compare handles the min internally
stats = compare(m, p_star, result.demands[:,1], result.supplies[:,1])
# stats.zi_efficiency uses min.(zi_d, zi_s) — actual trades, not gross demand
```

### Phase 7 — Padé Approximation for Irrational Utility Classes

Load via `include(joinpath(pkgdir(DiscreteMarket), "src", "extensions", "pade.jl"))`.

**Cobb-Douglas** demand is exactly rational: `xⱼ = floor(αⱼ · I / pⱼ)` where `I = p · ω`.

**CES demand** `xⱼ = I · αⱼ^σ · pⱼ^(−σ) / Σₗ αₗ^σ · pₗ^(1−σ)` involves `pⱼ^(−σ)`, which is irrational for non-integer σ. It is approximated by a [m/m] diagonal Padé approximant: a rational function `P(t)/Q(t)` in the scaled variable `t = (p − p₀)/p₀`, where `p₀` is the domain midpoint and coefficients are rationalized to denominator bound `Q_bound`.

```julia
include(joinpath(pkgdir(DiscreteMarket), "src", "extensions", "pade.jl"))

# Cobb-Douglas: exact rational arithmetic
cd = CobbDouglas([1//2, 1//2])
demand(cd, Price[4//1, 6//1], Int[10, 8])       # => [11, 7]
walras_residual(cd, Price[4//1, 6//1], Int[10, 8])  # => 2//1  (unspent budget)

# CES with σ=1.5 (irrational): Padé order 2, domain [1,8]
ces = CESApprox([1//2, 1//2], 1.5, 1//1, 8//1; m=2, Q_bound=1000)
demand(ces, Price[2//1, 5//1], Int[15, 10])     # => [24, 6]
```

Approximation accuracy for the canonical `p^(−σ)` function:

| Approximant | Domain | Order | Max relative error |
|---|---|---|---|
| `p^(−2)` (rational) | [1, 9] | 2 | < 1 × 10⁻¹⁰ (exact) |
| `p^(−0.5)` | [1, 9] | 2 | 1.6% |
| `p^(−1.5)` | [1, 4] | 3 | 0.32% |

The identity `p^(1−σ) = p · p^(−σ)` is used in the CES denominator to avoid constructing a second Padé (which would degenerate for integer 1−σ).

### Phase 8 — Stern-Brocot Focal Price Explorer

Load via `include(joinpath(pkgdir(DiscreteMarket), "src", "extensions", "stern_brocot.jl"))`.

The **Stern-Brocot depth** of a rational price `p/q` is the number of left/right steps to reach it in the Stern-Brocot tree, equal to `sum(continued_fraction(p/q)) − 1`. Low-depth prices are "arithmetically simple" — they arise naturally as focal points.

```julia
include(joinpath(pkgdir(DiscreteMarket), "src", "extensions", "stern_brocot.jl"))

sb_depth(1//1)   # => 0   (root)
sb_depth(1//2)   # => 1
sb_depth(2//3)   # => 2
sb_depth(3//5)   # => 3
sb_depth(7//5)   # => 4   (cf = [1,2,2])
```

**`simplest_rational(lo, hi)`** returns the shallowest Stern-Brocot node in `[lo, hi]` — the rational with the smallest denominator. This runs in O(depth) via a recursive interval-halving algorithm.

**`focal_price(m, ε)`** finds the simplest price that approximately clears the market:

```julia
# Market with p* = 7//5 (depth 4).  With ε = 1//5, the interval [6/5, 8/5]
# contains 3//2 (depth 2), which is the focal price.
focal_price(m, 1//5)   # => 3//2   (depth 2, vs WGE depth 4)

stats = focal_stats(m, 1//5)
# stats.depth_wge    => 4
# stats.depth_focal  => 2
# stats.depth_saving => 2
```

Over 200 random markets (2–6 consumers, 2–6 firms, up to 4 units, Q=100):

| Depth metric | Value |
|---|---|
| Minimum WGE depth | 0 |
| Mean WGE depth | 11.3 |
| Maximum WGE depth | 49 |

The depth_saving is always ≥ 0: the focal price is never more complex than the WGE price.

---

## Running the Tests

```
julia --project=. -e 'import Pkg; Pkg.test()'
```

The suite runs 37 test sets including three randomized batches of 200 markets each and phase-gated extension tests (phases 7–8). Four non-obvious findings surfaced during development:

**Equilibrium non-existence in multi-unit markets.** The existence theorem (`max_i vᵢ¹ ≥ min_f cᶠ¹ → equilibrium exists`) holds only for single-unit agents. With multi-unit demands, a demand step can jump by 2 (two consumers share the same WTP), skipping over the supply level entirely. `find_equilibrium` correctly returns `nothing` in those cases — 1 out of 200 random markets has no competitive equilibrium.

**Overflow in naïve rational bisection.** The midpoint `(p_lo + p_hi) / 2` doubles the rational denominator at each iteration, overflowing `Int64` within ~60 steps from the default bounds `[1/100, 1000]`. `tatonnement` uses the Stern-Brocot mediant `(a+c)//(b+d)` instead, which keeps denominators growing linearly.

**Wrong surplus formula.** The spec's `total_surplus` gave each consumer up to `q` units rather than distributing `q` units efficiently across all agents — overcounting consumer surplus by a factor proportional to the number of consumers. The correct formula sorts all WTP descending and all WTA ascending, then differences the top-`q` prefixes.

**`compare` used gross demand for efficiency.** Passing raw ZI demand (rather than `min(demand, supply)`) to `zi_efficiency` assumed every demanded unit was fulfilled. `compare` now applies `min.(zi_d, zi_s)` internally, and `generate_market` propagates `nothing` equilibria rather than silently substituting `1//1`.

---

## Roadmap

| Phase | Status | Description |
|-------|--------|-------------|
| 1–6 | Done | Core types, solvers, ZI simulation, comparison stats, multi-good |
| 7 | Done (extension) | Padé approximation for CES and irrational utility classes |
| 8 | Done (extension) | Stern-Brocot focal price explorer |
| 9 | Stub | Factor markets with endogenous wage `w` |
| 10 | Pending | OLG dynamic wrapper |
