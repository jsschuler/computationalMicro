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

Phases 1–5 of the [design spec](discrete_market_spec.md) are complete and tested.

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

mkt, p_stars = generate_market(42; k=3, n_consumers=4, n_firms=3, max_units=5)
```

### Phase 3 — WGE Solvers

**Exact scan** (`find_equilibrium`) — scans the sorted union of all WTP/WTA values. `O(n log n)`.

**Tatônnement** (`tatonnement`) — bisection on ℚ using the Stern-Brocot mediant to keep denominators bounded. Returns an exact rational equilibrium price.

```julia
find_equilibrium(m)    # => 5//1  (scan)
tatonnement(m)         # => 5//1  (bisection)
solve_wge_exact(mkt)   # => vector of prices, one per good
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

```julia
total_surplus(m, 5)            # => 20//1  (exact rational, efficient allocation)
wge_surplus(m, 5//1)           # => 20//1

zi_efficiency(m, 5//1, traded) # ∈ [0, 1]

stats = compare(m, p_star, zi_demands, zi_supplies)
# ComparisonStats: good, p_wge, q_wge, zi_mean_demand, zi_mean_supply,
#                 zi_efficiency, n_trials
```

---

## Running the Tests

```
julia --project=. -e 'import Pkg; Pkg.test()'
```

The suite runs 12 test sets including three randomized batches of 200 markets each. Two non-obvious findings surfaced:

**Equilibrium non-existence in multi-unit markets.** The existence theorem (`max_i vᵢ¹ ≥ min_f cᶠ¹ → equilibrium exists`) holds only for single-unit agents. With multi-unit demands, a demand step can jump by 2 (two consumers share the same WTP), skipping over the supply level. `find_equilibrium` correctly returns `nothing` in those cases — 1 out of 200 random markets has no competitive equilibrium.

**Overflow in naïve rational bisection.** The midpoint `(p_lo + p_hi) / 2` doubles the rational denominator at each iteration, overflowing `Int64` within ~60 steps from the default bounds `[1/100, 1000]`. `tatonnement` uses the Stern-Brocot mediant `(a+c)//(b+d)` instead, which keeps denominators growing linearly.

---

## Roadmap

| Phase | Status | Description |
|-------|--------|-------------|
| 1–5 | Done | Core types, solvers, ZI simulation, comparison stats |
| 6 | Pending | Full k-good `Market` with vector tatônnement |
| 7 | Pending | Padé approximation for CES and irrational utility classes |
| 8 | Pending | Stern-Brocot focal price explorer |
| 9 | Pending | Factor markets with endogenous wage `w` |
| 10 | Pending | OLG dynamic wrapper |
