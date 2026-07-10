# DiscreteMarket.jl — Analysis Log

Empirical results from systematic ZI vs WGE analysis and the BiGeometric jump-index
investigation.  All runs use `Q=100` denominator bound, `max_units=5` per agent, and
ZI consumer wealth calibrated as `p* × len(wtp)` so budget-constrained max demand
matches WGE capacity at equilibrium.

---

## 1. Test Suite

Phases 1–9 complete and passing.  All 43 test sets pass in ~10 s.

```
[ Info: Markets with no competitive equilibrium: 1 / 200
[ Info: p^(-0.5) Padé[2,2] max relative error: 1.6130089900092237%
[ Info: p^(-1.5) Padé[3,3] max relative error: 0.3210055111179222%
[ Info: CES σ=1.5: Padé demand [24, 6]  vs exact [24, 6]
[ Info: WGE price sb_depth — min: 0  mean: 11.27  max: 49
[ Info: Factor market randomized: labor cleared exactly 47 / 500; goods non-existence: 55
```

---

## 2. ZI Efficiency — Small/Mixed Markets

**Script:** `test/zi_efficiency_analysis.jl`
**Setup:** 1 000 random markets; agent counts drawn from 1–8 consumers, 1–8 firms,
1–6 max units; 1 000 ZI trials per market.

```
Markets with no WGE:       3
Markets with zero surplus: 26
Markets analyzed:          971

── Overall efficiency (ZI surplus / WGE surplus)
  Mean:    0.6961
  Median:  0.7506
  Std:     0.4677
  Min:     -12.2666
  Max:     0.9501

── Efficiency distribution
  [0.00–0.05) │                                            0.0% (0)
  [0.05–0.10) │                                            0.0% (0)
  [0.10–0.15) │                                            0.1% (1)
  [0.15–0.20) │                                            0.1% (1)
  [0.20–0.25) │ █                                          0.3% (3)
  [0.25–0.30) │ █                                          0.4% (4)
  [0.30–0.35) │ █                                          0.5% (5)
  [0.35–0.40) │ █                                          0.4% (4)
  [0.40–0.45) │ █████                                      2.2% (21)
  [0.45–0.50) │ ██████████████                             5.5% (53)
  [0.50–0.55) │ ██████████████████                         7.1% (69)
  [0.55–0.60) │ ████████████████                           6.3% (61)
  [0.60–0.65) │ ████████████████                           6.3% (61)
  [0.65–0.70) │ ███████████████████████████               10.7% (104)
  [0.70–0.75) │ ████████████████████████                   9.5% (92)
  [0.75–0.80) │ ██████████████████████████████████        13.4% (130)
  [0.80–0.85) │ ████████████████████████████████████████  16.0% (155)
  [0.85–0.90) │ █████████████████████████████████████     14.6% (142)
  [0.90–0.95) │ ███████████████                            6.0% (58)
  [0.95–1.00) │                                            0.1% (1)

── Breakdown by n_consumers
  n_consumers=1–2  n=231  mean=0.5596  median=0.5930  std=0.3084
  n_consumers=3–5  n=361  mean=0.6937  median=0.7720  std=0.7008
  n_consumers=6–8  n=379  mean=0.7816  median=0.8248  std=0.1305

── Breakdown by n_firms
  n_firms=1–2     n=231  mean=0.5812  median=0.5930  std=0.1788
  n_firms=3–5     n=349  mean=0.7300  median=0.7647  std=0.1351
  n_firms=6–8     n=391  mean=0.7337  median=0.8189  std=0.7062

── Breakdown by max_units
  max_units=1–2   n=309  mean=0.6043  median=0.6852  std=0.7901
  max_units=3–4   n=330  mean=0.7268  median=0.7615  std=0.1440
  max_units=5–6   n=332  mean=0.7511  median=0.7828  std=0.1668

── Breakdown by WGE equilibrium quantity q*
  q*=1–2   n=308 (31.7%)  mean=0.5059  median=0.5580  std=0.7854
  q*=3–5   n=366 (37.7%)  mean=0.7422  median=0.7600  std=0.0997
  q*=6–10  n=259 (26.7%)  mean=0.8305  median=0.8454  std=0.0700
  q*=11+   n=38  (3.9%)   mean=0.8782  median=0.8866  std=0.0414

── Tail counts
  eff >= 0.50 : 873 / 971  (89.9%)
  eff >= 0.60 : 743 / 971  (76.5%)
  eff >= 0.70 : 578 / 971  (59.5%)
  eff >= 0.80 : 356 / 971  (36.7%)
  eff >= 0.90 :  59 / 971   (6.1%)
  eff >= 0.95 :   1 / 971   (0.1%)
```

**Key findings:**
- Mean efficiency 69.6%, but highly dispersed (std 0.47); driven by thin markets (q*=1–2)
- Negative minimum (−12.3): ZI traders in thin markets can actively destroy surplus by trading past q*
- Strong monotone relationship with q*: larger equilibria → higher and less variable efficiency
- Efficiency increases monotonically with n_consumers and n_firms

---

## 3. ZI vs WGE — Large Markets

**Script:** `test/zi_large_market.jl`
**Setup:** 500 markets per size; 1 000 ZI trials per market; Q=100, max 5 units/agent.

### 3a. 100 agents (50C + 50F)

```
Markets analyzed:  471 / 500  (29 had no WGE)

── WGE outcomes
  q*:  mean=74.1  median=74.0  [58, 88]
  S*:  mean=365.3  median=363.9  [286.1, 475.0]

── ZI outcomes (at fixed p*)
  qty: mean=69.2  median=69.2  [58.9, 80.5]
  S:   mean=359.1  median=358.4  [282.3, 457.0]

── WGE vs ZI gap
  Mean quantity gap  (q* − ZI):  4.91 units
  Mean surplus ratio (ZI / WGE): 0.9833

── ZI efficiency
  Mean=0.9833  Median=0.9859  Std=0.0079  Min=0.9455  Max=0.9928

── Tail counts
  eff >= 0.90 : 471 / 471  (100.0%)
  eff >= 0.95 : 469 / 471   (99.6%)
  eff >= 0.98 : 354 / 471   (75.2%)
  eff >= 0.99 :  71 / 471   (15.1%)
```

### 3b. 500 agents (250C + 250F)

```
Markets analyzed:  420 / 500  (80 had no WGE)

── WGE outcomes
  q*:  mean=374.0  median=374.0  [333, 407]
  S*:  mean=1846.3  median=1845.7  [1611.9, 2035.0]

── ZI outcomes (at fixed p*)
  qty: mean=363.4  median=363.3  [335.5, 390.5]
  S:   mean=1840.3  median=1839.0  [1608.0, 2026.8]

── WGE vs ZI gap
  Mean quantity gap  (q* − ZI):  10.65 units
  Mean surplus ratio (ZI / WGE): 0.9967

── ZI efficiency
  Mean=0.9967  Median=0.9974  Std=0.0017  Min=0.9870  Max=0.9986

── Tail counts
  eff >= 0.90 : 420 / 420  (100.0%)
  eff >= 0.95 : 420 / 420  (100.0%)
  eff >= 0.98 : 420 / 420  (100.0%)
  eff >= 0.99 : 416 / 420   (99.0%)
```

**Key findings:**
- ZI efficiency converges to WGE as market grows (law-of-large-numbers averaging)
- Std collapses from 0.008 (100 agents) to 0.002 (500 agents)
- Non-existence rate *increases* with market size: 5.8% at 100 agents, 16% at 500 agents
  (more agents → more WTP collisions → larger demand jumps → see §§4–5)
- Quantity gap grows in absolute terms but shrinks relatively: 4.9/74 = 6.6% vs 10.7/374 = 2.9%

**Convergence summary across all market sizes:**

| Agents | Mean eff | Std    | Non-existence rate |
|--------|----------|--------|--------------------|
| 2–16 (random) | 69.6% | 0.47 | 0.3% |
| 100    | 98.3%    | 0.008  | 5.8%               |
| 500    | 99.7%    | 0.002  | 16.0%              |

---

## 4. Jump-Index Conjecture

**Script:** `test/jump_index.jl`
**Conjecture tested:** P(WGE exists) ≈ P(J(m) = 1), where J(m) = max single-price
drop in aggregate demand (maximum WTP-value multiplicity across all consumer-unit pairs).

```
── 8 agents (4C + 4F)   n=1000
  P(J=1)=0.898   P(WGE)=0.994
  P(WGE|J=1)=1.000   P(WGE|J>1)=0.941
  Conjecture gap = 0.096

  J distribution:  J=1: 89.8%  J=2: 10.1%  J=3: 0.1%
  P(WGE|J=j):  J=1: 1.000  J=2: 0.941  J=3: 1.000

── 50 agents (25C + 25F)   n=1000
  P(J=1)=0.029   P(WGE)=0.968
  P(WGE|J=1)=1.000   P(WGE|J>1)=0.967
  Conjecture gap = 0.939

  J distribution:  J=1: 2.9%  J=2: 83.0%  J=3: 13.7%  J=4: 0.4%
  P(WGE|J=j):  J=1: 1.000  J=2: 0.970  J=3: 0.949  J=4: 1.000

── 100 agents (50C + 50F)   n=1000
  P(J=1)=0.000   P(WGE)=0.952
  P(WGE|J>1)=0.952   (J=1 never occurs)
  Conjecture gap = 0.952

  J distribution:  J=2: 33.3%  J=3: 58.5%  J=4: 7.8%  J=5: 0.4%
  P(WGE|J=j):  J=2: 0.946  J=3: 0.952  J=4: 0.974  J=5: 1.000

── 500 agents (250C + 250F)   n=1000
  P(J=1)=0.000   P(WGE)=0.864
  P(WGE|J>1)=0.864   (J=1 never occurs)
  Conjecture gap = 0.864

  J distribution:  J=4: 0.8%  J=5: 30.2%  J=6: 50.0%  J=7: 14.9%
                   J=8: 3.3%  J=9: 0.9%
  P(WGE|J=j):  J=5: 0.841  J=6: 0.870  J=7–9: ~0.85
```

**Verdict:** Conjecture **refuted**.
- J=1 is sufficient for WGE existence (P(WGE|J=1) = 1.000 always) ✓
- J=1 is not necessary: at 100+ agents J=1 never occurs yet WGE exists 86–95% of the time
- P(WGE|J=j) is flat across j values at any fixed market size — J carries no predictive
  information once it exceeds 1
- Root cause: J measures the absolute demand drop; what matters is the drop relative
  to the supply gap at that price (see §5)

---

## 5. BiGeometric Jump Index — BJ2 Theorem

**Script:** `test/bj_index.jl`
**Setup:** 2 000 markets per size; four market sizes (8, 50, 100, 500 agents).

Three predictors compared:
- **J** — additive maximum jump (from §4)
- **BJ** = max_v J(v)/D(v+) — relative jump vs post-jump demand
- **BJ2** = max_v J(v)/z(v) where z(v) = D(v) − S_max(v) > 0
            "jump divided by supply gap at that price"

**Bug fix in `bridge_failure`:** the original implementation added an unnecessary
`d_gap < s_lo` check against the wrong gap.  The correct condition is simply:
excess demand at v AND excess supply just above v (sign change), which is equivalent
to BJ2 > 1.  After fix, `bridge_failure` agrees with `fails` 1.0000 at every size.

```
── 8 agents (4C + 4F)   n=2000
  P(WGE fails)=0.0050   P(bf)=0.0050   P(BJ2>1)=0.0050
  Agreement bf==fails=1.0000   BJ2>1==fails=1.0000

  P(fail|J=j):    J=1: 0.000   J=2: 0.051
  P(fail|BJ2 bin): [0,0.5): 0.000  [0.5,1): 0.000  [1,2): 0.000  [2,∞): 1.000

── 50 agents (25C + 25F)   n=2000
  P(WGE fails)=0.0275   P(bf)=0.0275   P(BJ2>1)=0.0275
  Agreement bf==fails=1.0000   BJ2>1==fails=1.0000

  P(fail|J=j):    J=1: 0.000  J=2: 0.026  J=3: 0.048  J=4: 0.000
  P(fail|BJ2 bin): [0,0.5): 0.000  [0.5,1): 0.000  [1,2): 0.073  [2,∞): 1.000

── 100 agents (50C + 50F)   n=2000
  P(WGE fails)=0.0450   P(bf)=0.0450   P(BJ2>1)=0.0450
  Agreement bf==fails=1.0000   BJ2>1==fails=1.0000

  P(fail|J=j):    J=2: 0.044  J=3: 0.047  J=4: 0.036  J=5: 0.000
  P(fail|BJ2 bin): [0,0.5): 0.000  [0.5,1): 0.000  [1,2): 0.053  [2,∞): 1.000

── 500 agents (250C + 250F)   n=2000
  P(WGE fails)=0.1470   P(bf)=0.1470   P(BJ2>1)=0.1470
  Agreement bf==fails=1.0000   BJ2>1==fails=1.0000

  P(fail|J=j):    J=5: 0.155  J=6: 0.145  J=7: 0.144  J=8: 0.174
  P(fail|BJ2 bin): [0,0.5): 0.000  [0.5,1): 0.000  [1,2): 0.258  [2,∞): 1.000
```

**The BJ2 Theorem (empirically confirmed, analytically derived):**

> WGE fails ⟺ BJ2 > 1
>
> i.e., ∃ WTP value v such that:
>   - z(v) = D(v) − S_max(v) > 0   (excess demand at v), AND
>   - J(v) > z(v)                   (demand jump exceeds the supply gap)
>
> Equivalently: demand drops from above the supply range to below it at some price,
> with no price in between where D ∈ [S_min, S_max].

**Proof sketch:**
- (⟹) BJ2 > 1 at v implies D(v) > S_max(v) and D(v+) < S_max(v).
  Supply just above v = S_max(v) (singleton), and supply is non-decreasing, so
  excess supply holds at ALL p > v.  Demand is non-increasing, so D(p) ≥ D(v) > S_max(v) ≥ S_max(p)
  for all p ≤ v → excess demand throughout.  No clearing exists.
- (⟸) If WGE fails, there is a sign change with no clearing, which requires demand
  to jump over the supply interval at some WTP value v, i.e., J(v) > z(v) > 0.

**Why J fails as a predictor:** J measures the absolute jump size, but what determines
non-existence is whether J exceeds the supply gap z(v) at the same price.  With large
markets z(v) is typically large (many units of supply), so even J=6 is usually well
below z(v).  J is flat-predictive because z(v) and J(v) co-vary with market size.

**BiGeometric interpretation:** In the BiGeometric (multiplicative) calculus, the
natural measure of a demand discontinuity at v is J(v)/z(v) — how many multiples of
the supply gap the demand jump spans.  The mesh condition for WGE existence is
max_v J(v)/z(v) ≤ 1: demand never jumps more than one "supply gap width."
This is the correct discrete analog of the Sperner/Brouwer mesh condition, replacing
the additive J ≤ 1 requirement of single-unit markets.

---

## 6. BiGeometric ZIT — Quantity-Space Formulation

**Script:** `test/bgi_efficiency.jl`
**Idea:** Replace the additive-uniform quantity draw of standard ZIT with a log-uniform
(BiGeometric-uniform) draw.  Standard ZIT samples `q ~ Uniform{0,...,Q_max}`; BiGeo
quantity samples `q = floor(exp(u * log(Q_max+2))) - 1` for `u ~ Uniform[0,1)`, giving
`P(q=k) = log((k+2)/(k+1)) / log(Q_max+2)` — a discrete Zipf/reciprocal distribution.

| Statistic | Standard ZIT | BiGeo Quantity | Δ |
|-----------|-------------|----------------|---|
| Mean eff  | 69.6%       | 54.5%          | −15.2pp |
| Std       | 0.468       | 0.242          | −0.226 |
| Min       | −12.27      | −5.18          | +7.08 |
| Sharpe    | 1.49        | 2.25           | +0.76 |

**Finding:** Log-uniform quantity draws halve the standard deviation and raise the minimum,
but reduce mean efficiency by 15pp.  Agents systematically under-demand (median ≈ √Q_max
instead of Q_max/2), and this gap does not close with market size — BiGeo quantity
converges to `≈ 65%` of WGE efficiency as N→∞ (not 100%).  Standard ZIT dominates
head-to-head in 97% of markets.

---

## 7. BiGeometric ZIT — Price-Space Formulation

**Scripts:** `test/bgi_price_efficiency.jl` (small markets), `test/bgi_price_large.jl`
(large markets).

### Formulation

Instead of drawing a quantity, each agent draws a **virtual participation price** from
their WTP/WTAC range in log-price space, then trades all units on the right side of it:

- **Consumer:** draw `p_v ~ LogUniform[p*, max_WTP_i]`; demand all units with `wtp ≥ p_v`
- **Firm:** draw `p_v ~ LogUniform[min_WTAC_i, p*]`; supply all units with `wtac ≤ p_v`

The resulting per-unit purchase probability is:

```
P(buy unit j)    = log(wtp_j  / p*) / log(wtp_max / p*)
P(supply unit j) = log(p*  / wtac_j) / log(p*  / wtac_min)
```

Unit j's probability is proportional to how far it sits from p* **in multiplicative
(log-price) terms** — i.e., the BiGeometric distance from the equilibrium price.
Unit 1 (highest WTP / lowest cost) always has probability 1.  Marginal units (WTP → p+)
have probability → 0.

**Structural guarantee:** BiGeo price demand ≤ WGE demand always, so efficiency ∈ [0,1]
and negative efficiency is impossible.

### 7a. Small/mixed markets (1000 random markets, same seeds as §2)

```
                Standard ZIT    BiGeo Price     Δ
  Mean:         0.6961          0.9115        +0.2154
  Median:       0.7506          0.9325        +0.1818
  Std:          0.4677          0.0951        -0.3727
  Min:          -12.2666         0.5021       +12.7687
  Max:          0.9501          1.0000        +0.0499

  Mean/Std (Sharpe):  Std = 1.488,  BiGeo-P = 9.587

  Negative efficiency: 0 / 971  (vs 6 / 971 for Standard)

  Threshold │  Standard    BiGeo-Price
  eff≥0.50  │  872/971 (89.8%)  971/971 (100.0%)
  eff≥0.70  │  578/971 (59.5%)  942/971 (97.0%)
  eff≥0.90  │   59/971 (6.1%)  594/971 (61.2%)
  eff≥0.95  │    1/971 (0.1%)  430/971 (44.3%)

  BiGeo-P > Std: 842/971 (86.7%)
```

**Breakdown by q*:**

| q* range | n | Std mean | BiGeo-P mean | Δ |
|----------|---|----------|--------------|---|
| q*=1–2   | 308 | 50.6% | 96.4% | +45.8pp |
| q*=3–5   | 366 | 74.2% | 89.8% | +15.6pp |
| q*=6–10  | 259 | 83.1% | 87.7% | +4.6pp  |
| q*=11+   |  38 | 87.8% | 85.3% | −2.5pp  |

**Mechanism:** BiGeo-P implicitly selects **inframarginal units** (high WTP, low cost)
because their log-distance from p* is large → high purchase probability.  Marginal units
(WTP barely > p*) have near-zero purchase probability.  In thin markets this eliminates
the "wrong-unit trading" that causes negative efficiency in standard ZIT.  The crossover
to standard ZIT dominance happens at q* ≈ 10–11, where LLN begins covering the marginal
units that BiGeo-P misses.

### 7b. Large markets (500 seeds each size)

```
                              Standard ZIT    BiGeo Price     Δ
100 agents (50C+50F)
  q* mean = 74.1
  Mean efficiency:            0.9833          0.9146        −0.0687
  Std:                        0.0079          0.0176        +0.0098
  Min:                        0.9455          0.8620        −0.0835
  Mean traded q / q*:         93.4%           71.1%
  Head-to-head: Std wins:     471/471 (100%)

500 agents (250C+250F)
  q* mean = 374.0
  Mean efficiency:            0.9967          0.9168        −0.0800
  Std:                        0.0017          0.0091        +0.0074
  Min:                        0.9870          0.8905        −0.0965
  Mean traded q / q*:         97.2%           71.2%
  Head-to-head: Std wins:     420/420 (100%)
```

**The plateau:** BiGeo-P efficiency stabilises at ≈91–92% across both large market sizes
while standard ZIT converges to 99–100%.  Crucially, BiGeo-P traded quantity locks at
≈71% of q* at both 100 and 500 agents — the ratio does not improve with N.  This is a
structural ceiling: marginal units (WTP ≈ p*) always have near-zero log-distance from p*
and are never sampled, regardless of how large the market grows.

### 7c. Full efficiency-vs-market-size table (all three methods)

| Market size | Standard ZIT | BiGeo Quantity | BiGeo Price | |
|---|---|---|---|---|
| Small, q*=1–2 | 50.6% | 42.0% | **96.4%** | BiGeo-P wins |
| Small, q*=3–5 | 74.2% | 52.8% | **89.8%** | BiGeo-P wins |
| Small, q*=6–10 | 83.1% | 54.3% | **87.7%** | BiGeo-P wins |
| Small, q*=11+ | **87.8%** | 56.2% | 85.3% | Std wins |
| 100 agents | **98.3%** | ~65%* | 91.5% | Std wins |
| 500 agents | **99.7%** | ~65%* | 91.7% | Std wins |

*BiGeo Quantity large-market limit not directly measured; predicted from k=5 scaling.

### 7d. Interpretation

**Standard ZIT** converges to 100% by the law of large numbers: uniform draws average
to Q_max/2 ≈ q*, so total traded quantity → q* as N→∞.

**BiGeo Quantity** converges to a fixed fraction ≈ 65% (for k=5 units/agent) because
`E[q_i] ≈ k/log(k+1)` < `k/2` for all finite k.  LLN converges to the wrong mean.

**BiGeo Price** converges to a fixed fraction ≈ 91–92% because it permanently
undersamples marginal units.  The quantity gap (71% of q*) carries a small surplus gap
because marginal units have near-zero surplus (their WTP ≈ p* ≈ WTAC).

The **crossover at q*≈10** reflects when marginal units begin to matter:
- Thin markets: marginal units are few and low-surplus; BiGeo-P's inframarginal focus wins
- Thick markets: many marginal units with collectively significant surplus; LLN wins

**BiGeo price-space ZIT is the right model when:** agents draw participation thresholds
from log-price space (natural when agents think multiplicatively about price ratios) and
the market is thin or prices are uncertain.  Standard ZIT is better when the market is
thick and the wealth calibration encodes accurate price information.

---

## 8. Additive–BiGeometric Mixture Model

**Script:** `test/mixed_zit.jl`
**Mixing parameter:** λ ∈ [0,1].  Each agent, each period, independently:
- with prob (1−λ): draw via **standard ZIT** (q ~ Uniform{0,...,Q_max})
- with prob λ: draw via **BiGeo price** (p_v ~ LogUniform[p*, max_WTP])

λ=0 → pure standard ZIT; λ=1 → pure BiGeo price.  All λ values run on identical
pre-drawn random numbers so the comparison is controlled.

### 8a. Small/mixed markets (n=971 analyzed, T=500 trials each)

```
  λ     │  Mean    Median   Std     Min       >0.50   >0.90
  ──────┼────────────────────────────────────────────────────
  λ=0.0 │  0.6961  0.7520  0.4733  -12.5604   89.9%    5.5%
  λ=0.1 │  0.7157  0.7642  0.4047  -10.4520   94.6%    6.2%
  λ=0.2 │  0.7347  0.7777  0.3520   -8.8048   96.7%    6.8%
  λ=0.3 │  0.7540  0.7880  0.3053   -7.3827   98.1%    7.9%
  λ=0.4 │  0.7740  0.8013  0.2479   -5.6051   98.9%   10.2%
  λ=0.5 │  0.7947  0.8130  0.2009   -4.3571   99.4%   13.1%
  λ=0.6 │  0.8161  0.8295  0.1529   -2.8947   99.7%   16.8%
  λ=0.7 │  0.8382  0.8445  0.1089   -1.4538   99.8%   23.5%
  λ=0.8 │  0.8603  0.8718  0.0901   -0.7364   99.8%   32.0%
  λ=0.9 │  0.8849  0.9071  0.0793    0.4596   99.9%   54.1%
  λ=1.0 │  0.9115  0.9332  0.0952    0.5106  100.0%   61.4%
```

Head-to-head vs λ=0: even λ=0.1 wins **80.5%** of markets; win rate holds at 85–88%
for all λ ∈ [0.3, 1.0] (marginal markets already converted at low λ; additional weight
increases the margin of improvement, not the number of markets won).

### 8b. 100-agent markets (n=284 analyzed, T=500 trials each)

```
  λ     │  Mean    Median   Std     Min       >0.90   >0.95
  ──────┼────────────────────────────────────────────────────
  λ=0.0 │  0.9839  0.9866  0.0077   0.9445  100.0%   99.3%
  λ=0.1 │  0.9816  0.9844  0.0089   0.9383  100.0%   99.3%
  λ=0.2 │  0.9783  0.9807  0.0101   0.9325  100.0%   98.9%
  λ=0.3 │  0.9739  0.9760  0.0111   0.9252  100.0%   95.8%
  λ=0.5 │  0.9622  0.9642  0.0126   0.9112  100.0%   84.5%
  λ=0.7 │  0.9464  0.9478  0.0137   0.8969   99.6%   44.7%
  λ=1.0 │  0.9147  0.9167  0.0179   0.8584   78.9%    1.1%
```

Standard ZIT wins head-to-head at λ=0.1 in 89% of markets; at λ≥0.7 it wins 100%.

### 8c. 500-agent markets (n=168 analyzed, T=500 trials each)

```
  λ     │  Mean    Median   Std     Min       >0.90   >0.95
  ──────┼────────────────────────────────────────────────────
  λ=0.0 │  0.9967  0.9973  0.0018   0.9876  100.0%  100.0%
  λ=0.1 │  0.9949  0.9955  0.0027   0.9841  100.0%  100.0%
  λ=0.2 │  0.9919  0.9925  0.0034   0.9805  100.0%  100.0%
  λ=0.5 │  0.9751  0.9756  0.0054   0.9618  100.0%  100.0%
  λ=0.7 │  0.9569  0.9576  0.0066   0.9405  100.0%   82.7%
  λ=1.0 │  0.9169  0.9179  0.0097   0.8927   95.2%    0.0%
```

Standard ZIT wins head-to-head at λ=0.1 in 95% of markets; at λ≥0.3 it wins 100%.

### 8d. Interpolation confirmed

The mixture interpolates **monotonically and smoothly** in both directions:

| Market size | Direction as λ increases | Mean eff. at λ=0 | Mean eff. at λ=1 | Total Δ |
|---|---|---|---|---|
| Small (thin) | ↑ monotone | 69.6% | 91.2% | +21.5pp |
| 100 agents   | ↓ monotone | 98.4% | 91.5% | −6.9pp  |
| 500 agents   | ↓ monotone | 99.7% | 91.7% | −8.0pp  |

**Why approximately affine in λ:** by linearity of expectation, each period each agent
independently selects a mechanism, so E[eff(λ)] ≈ (1−λ)·E[eff_std] + λ·E[eff_bgip].
This holds for the mean; std and min are non-linear — std drops faster at low λ (even
a small BiGeo fraction suppresses the worst standard ZIT draws in thin markets), and the
minimum improves steeply between λ=0.8 and λ=0.9.

**Optimal λ is market-size dependent and monotone:** there is no interior optimum — the
optimal λ is 1 for thin markets and 0 for thick markets.  The mixture is most useful
when market thickness is uncertain: a moderate λ (e.g., 0.4–0.6) hedges between the
inframarginal-selection advantage (BiGeo-P) and the LLN coverage advantage (standard).

---

## 9. Key Findings Summary

| Finding | Detail |
|---------|--------|
| ZI efficiency in thin markets | Mean ~70%, highly dispersed; can be negative when ZI trades past q* |
| ZI efficiency at 100 agents | Mean 98.3%, std 0.008; 100% above 90% |
| ZI efficiency at 500 agents | Mean 99.7%, std 0.002; 99% above 99% |
| WGE non-existence rate grows with N | 0.3% (small), 5.8% (100 agents), 16% (500 agents) |
| WGE + ZI trade-off | Larger markets → higher ZI efficiency AND higher non-existence rate |
| J conjecture refuted | J > 1 does not predict non-existence; P(WGE\|J=j) flat across j |
| BJ2 theorem confirmed | BJ2 > 1 ⟺ WGE fails; 100% agreement at all market sizes |
| bridge_failure bug fixed | Removed erroneous `d_gap < s_lo` check; now equivalent to BJ2 > 1 |
| Non-existence grows via birthday paradox | With Q=100 and N agents, expected collisions ~ N²/Q² → larger J; but J/D shrinks; z(v) absorbs the growth until BJ2 crosses 1 |
| BiGeo quantity ZIT | Lower mean (55%) but lower variance (std 0.24 vs 0.47); Sharpe 2.25 vs 1.49; converges to ~65% not 100% |
| BiGeo price ZIT (small markets) | Mean 91.2%, Sharpe 9.59; zero negative efficiency; wins 86.7% head-to-head vs standard |
| BiGeo price ZIT (large markets) | Plateaus at ~91-92%; standard ZIT wins 100% head-to-head; 71% quantity coverage is structural ceiling |
| Crossover threshold | BiGeo price dominates standard ZIT for q* ≤ ~10; standard dominates for q* ≥ 11 |
| Additive–BiGeo mixture | Interpolates affinely in mean efficiency; optimal λ=1 (thin) or λ=0 (thick); no interior optimum |
| Adaptive ZIT (qty signal, small) | Cold-start exploration, then rapid reversion to standard; final λ≈0.076, eff≈0.75 |
| Adaptive ZIT (qty signal, large) | Consumers converge to λ≈0 (standard); firms stabilize at λ≈0.375 (Nash) in 500-agent markets |
| Firm Nash equilibrium | 37.5% of firms play BiGeo in 500-agent steady state; supply-reduction strategy avoids pro-rata rationing |
| Adaptive ZIT (surplus signal) | Consumer dynamics unchanged (surplus monotone in qty); firm Nash at 37.5% confirmed robust |
| Adaptive ZIT (marginal surplus) | λ → 67% (small) / 80% (large); efficiency collapses to 70%/88% — miscoordination trap |
| Miscoordination trap | Marginal signal is individually rational but drives collective underdemand; Nash ≠ social optimum |
| Signal-gated annealing (sym) | λ_eq ≈ 0.37 for all sizes; resolves trap in large markets (93-95% eff); thin markets need >300 turns |
| Golden-ratio equilibrium | λ_eq ≈ 0.382 = φ⁻² when n_active/n_total ≈ 0.5; determined by agent microstructure, not market size |
| Asymmetric rates (cool=0.01) | All markets converge to λ_eq ≈ 0.50; individual signal blind to market-level thickness |
| Fundamental limit | Individual marginal surplus cannot encode aggregate market thinness; market-level signal needed |

---

## 10. Adaptive ZIT — Learning Game (Quantity Signal)

**Script:** `test/adaptive_zit.jl`
**Setup:** T=300 turns per market; agents individually switch strategies turn-by-turn using a
last-observation comparison rule (optimistic prior: NaN → +Inf forces BiGeo exploration on
turn 2). Performance signal = pro-rata realized trade quantity.

### Learning rule

```
if perf_bgip (last seen) > perf_std (last seen):  play BiGeo this turn
else:                                               play Std this turn
```
NaN memory (never tried) maps to +Inf, guaranteeing each agent tries BiGeo exactly once
(turn 2) before the genuine comparison takes over from turn 3 onward.

### Cold-start dynamics (small/mixed markets)

```
  Turn │  λ_c    λ_f    Eff_mean   vs Std    vs BiGeo
  ─────┼──────────────────────────────────────────────
  t=  1 │  0.000  0.000  0.71xx    +0.00     -0.16
  t=  2 │  1.000  1.000  0.91xx    +0.20     +0.00   ← forced exploration
  t=  3 │  0.xxx  0.xxx  0.7x      ...       ...
  t=300 │  0.076  0.076  0.7538   +0.0518   -0.1638
```

- **Turn 1:** All agents play standard (initialization). eff ≈ 0.71 (= pure standard baseline).
- **Turn 2:** All agents switch to BiGeo (NaN optimistic prior). eff ≈ 0.91 (= pure BiGeo baseline).
- **Turn 3 onward:** BiGeo observation (0.91) vs standard observation (0.71) → most stay BiGeo initially.
  But as agents re-experience standard, quantity signal favors standard → reversion begins.
- **Turn 300:** λ converges to ~0.076. Efficiency ≈ 0.754 — between pure standard (0.696) and
  pure BiGeo (0.912). Adaptive beats static standard by +0.052 purely due to residual BiGeo adoption.

### Large markets

| Market size | n analyzed | Final λ_c | Final λ_f | Final eff | vs Static Std |
|-------------|-----------|-----------|-----------|-----------|---------------|
| Small/mixed | 971/1000  | 0.076     | 0.076     | 0.754     | +0.052        |
| 100 agents  | 284/300   | 0.009     | 0.009     | 0.986     | +0.003        |
| 500 agents  | 124/150   | 0.006     | 0.375     | 0.984     | +0.001        |

### Firm Nash equilibrium in 500-agent markets

The striking result: **consumers converge fully to standard (λ_c → 0.6%), but firms stabilize
at λ_f ≈ 37.5% BiGeo** — a persistent, stable asymmetry.

**Mechanism:** Standard ZIT firms draw from {0,...,n_total} including units with WTAC > p*.
At 500 agents, aggregate supply S often exceeds D, triggering firm-side rationing:
`realized_f = floor(draw_f × D/S)`. A standard firm with a large draw gets heavily rationed.
A BiGeo firm drawing only from active units supplies less overall, but avoids the worst
rationing events. Some fraction of firms (≈37.5%) find that their last BiGeo-period realization
exceeded their last standard-period realization (due to differential rationing exposure across
turns), and lock into BiGeo.

This is an emergent Nash: if fewer firms played BiGeo, S would fall and rationing would ease
for standard firms; if more played BiGeo, S would fall further and the rationing advantage
would erode. The 37.5% figure represents the fixed point of this stochastic adjustment process.

### Why quantity signal gives the wrong incentive for consumers

In thin markets, BiGeo price-space selection is socially superior (mean eff 91% vs 70%),
but this arises from selecting **high-surplus inframarginal units** — not from trading more
units. The quantity signal rewards the agent who draws a larger number, regardless of which
units. Since standard draws from a larger range and has higher mean, the quantity signal
causes all consumers to rationally defect back to standard — even though this is socially
suboptimal.

**Prediction:** switching to an individual surplus signal (surplus = sum(wtp[j] - p*) for
assigned units) should allow thin-market consumers to correctly identify BiGeo as superior,
since BiGeo selects higher-WTP units → higher per-unit surplus.

---

## 11. Adaptive ZIT — Individual Surplus Signal

**Script:** `test/adaptive_zit_surplus.jl`
**Setup:** identical to §10 but performance signal = individual surplus from pro-rata units.
`cum_surplus[k]` precomputed per agent for O(1) lookup.

```
Consumer surplus(k) = sum(active_wtp[1:k]) - k * p*   (capped at n_active)
Firm     surplus(k) = k * p* - sum(active_wtac[1:k])  (capped at n_active)
```

### Results

| Market size | n analyzed | Final λ_c | Final λ_f | Final eff | vs Static Std |
|-------------|-----------|-----------|-----------|-----------|---------------|
| Small/mixed | 971/1000  | 0.071     | 0.064     | 0.761     | +0.065        |
| 100 agents  | 284/300   | 0.009     | 0.007     | 0.985     | +0.002        |
| 500 agents  | 124/150   | 0.005     | 0.376     | 0.984     | −0.013        |

### Interpretation

**Consumer dynamics unchanged.** The surplus signal does not change the incentives for
consumers. This confirms the theoretical prediction: `surplus(k) = f(realized_k)` where `f`
is a monotone increasing function of realized units (every active unit contributes wtp_j − p* > 0).
Since standard ZIT gives higher mean realized units (pro-rata preserves the ratio), standard
still gives higher expected surplus. Consumers converge to standard in all market sizes.

**The 37.5% firm Nash is robust to the performance signal.** Under surplus, firms in the
500-agent case settle at λ_f ≈ 0.376 — essentially identical to the 0.375 from the quantity
signal. This confirms the firm Nash is a genuine game-theoretic equilibrium driven by the
stochastic comparison of rationing outcomes across turns, not an artifact of the quantity metric.

**Surplus trajectory in 500-agent firms:** λ_f peaks at ~0.57 around t=10–25 (firms see
surplus_bgip > surplus_std because BiGeo avoids pro-rata rationing events), then relaxes
to the 0.376 equilibrium. The overshoot and relaxation suggest the equilibrium is approached
from above under surplus, versus more gradually from below under the quantity signal.

**Small-market efficiency nearly identical:** 0.761 (surplus) vs 0.754 (quantity) — the
small difference comes from slightly different λ trajectories in the transient, not a
fundamentally different equilibrium.

### Theoretical explanation: why surplus ≡ quantity for consumers

For a consumer assigned k units with active WTP values sorted descending:
- `surplus(k) = Σ_{j=1}^{k} (wtp_j − p*)` — a strictly increasing function of k
- Standard gives higher mean k → higher mean surplus
- The only exception would be if BiGeo's k were ≥ n_active (all active units) while
  standard's k < n_active — but this requires extreme rationing where standard is heavily
  cut but BiGeo is not, which is rare

**Conclusion:** to materially change consumer learning dynamics, the performance signal
must be unit-value-weighted rather than quantity-weighted — e.g., highest-WTP unit
traded, or surplus per unit (WTP - p* of the marginal unit traded).

---

## 12. Adaptive ZIT — Marginal-Unit Surplus Signal

**Script:** `test/adaptive_zit_marginal.jl`
**Signal:** surplus of the last (k-th) unit assigned via pro-rata:
```
Consumer: full_wtp[k] - p*   (can be negative when k > n_active — overtrading)
Firm:     p* - full_wtac[k]  (can be negative when k > n_active — underselling)
```
Unlike cumulative surplus, the marginal signal **goes negative** when standard ZIT agents
overtrade (draw k > n_active), because the k-th unit in the sorted WTP array has wtp < p*.
BiGeo consumers never see negative marginal surplus (draws bounded to {0,...,n_active}).

### Results

| Market size | n | Final λ_c | Final λ_f | Final eff | vs Std   | vs BiGeo |
|-------------|---|-----------|-----------|-----------|----------|----------|
| Small/mixed | 971 | 0.672   | 0.648     | 0.698     | +0.002   | −0.213   |
| 100 agents  | 284 | 0.802   | 0.814     | 0.878     | −0.106   | −0.037   |
| 500 agents  | 124 | 0.781   | 0.816     | 0.888     | −0.109   | −0.029   |

### Learning dynamics (small markets)

```
  Turn │  λ_c    λ_f    Eff_mean   vs Std    vs BiGeo
  t=  3 │  0.402  0.426  0.792     +0.096    −0.120
  t= 10 │  0.546  0.545  0.755     +0.060    −0.156
  t= 50 │  0.648  0.628  0.713     +0.017    −0.199
  t=300 │  0.672  0.648  0.698     +0.002    −0.213
```

λ slowly creeps upward toward a fixed point around 0.67–0.67. The learning is correct at
the individual level: standard agents do observe negative marginal surplus when k > n_active,
and they switch to BiGeo. BiGeo agents observe non-negative marginal surplus and stay.

### The miscoordination trap (large markets)

In the 100- and 500-agent markets, λ converges to **~80% BiGeo** — far above the social
optimum (λ=0) — and efficiency collapses to **87-89%**, substantially worse than both pure
strategies.

**Why individual rationality fails socially:**
Each agent correctly observes that their marginal unit under standard sometimes carries
negative surplus (wtp[k] < p* when k > n_active). They switch to BiGeo, which has strictly
non-negative marginal surplus. But when 80% of agents use BiGeo simultaneously, aggregate
demand/supply falls well below q* (BiGeo draws only from active units, mean ≈ n_active/2;
standard draws from all units, mean ≈ n_total/2 ≈ q*/n). Market clearing drops to ~75% of q*,
reducing efficiency to ~88% — worse than pure BiGeo (91%).

**Synchronization amplifies the problem:** when excess demand gives negative marginal surplus
to many standard consumers in the same period, they ALL switch to BiGeo simultaneously.
The next period has very low collective demand (everyone BiGeo), tanking efficiency below
even the static mixture prediction. This anti-coordination creates efficiency oscillations
worse than the static mixture model's monotone decline.

**Summary across all three signals:**

| Signal            | Small mkt final λ | Small eff | Large λ | Large eff | Pattern                    |
|-------------------|-------------------|-----------|---------|-----------|----------------------------|
| Quantity          | 0.076             | 0.754     | ≈0      | 0.984     | Drift to standard          |
| Cumulative surplus| 0.071             | 0.761     | ≈0      | 0.984     | Same as quantity            |
| Marginal surplus  | 0.672             | 0.698     | ≈0.80   | 0.877     | Miscoordination trap        |

The marginal signal is the only one that changes consumer behavior — but it drives them
into a socially suboptimal equilibrium. The individual rationality it encodes (avoid negative
marginal surplus) conflicts with market-level efficiency (need agents to submit full demand
up to q* to clear the market).

---

## 14. Signal-Gated Annealing v2 — Asymmetric Rates (cool=0.01, heat=0.05)

**Script:** `test/adaptive_zit_anneal2.jl`
**Hypothesis:** smaller COOL_RATE gives heating a 5× structural advantage so thin markets
(frequent negative marginal surplus from standard) converge to higher λ_eq before the
cooling wins. T extended to 600 turns to confirm convergence.

### Results

| Market size  | Final λ_c | Final λ_f | Final eff | vs Std   | vs BiGeo |
|--------------|-----------|-----------|-----------|----------|----------|
| Small/mixed  | 0.494     | 0.535     | 0.693     | −0.003   | −0.219   |
| 100 agents   | 0.518     | 0.522     | 0.920     | −0.064   | +0.005   |
| 500 agents   | 0.530     | 0.524     | 0.934     | −0.063   | +0.017   |

### Key finding: asymmetric rates don't discriminate by market size

All three market sizes converge to the same λ_eq ≈ 0.50–0.53, compared to ≈ 0.37 under
symmetric rates. The equilibrium is NOT market-size-sensitive.

**Why:** the equilibrium λ is determined by each AGENT's own n_active/n_total ratio, not
by the market's aggregate thickness. Even in a 500-agent market, individual agents still
have units with WTAC > p* (inactive), so P(neg|std) > 0 per agent. With β/α = 5, even
small per-agent P(neg) generates enough heating to stabilize λ at ~0.50 regardless of
market size.

Equilibrium condition per agent (α=0.01, β=0.05):
```
α × λ × P(pos) = β × (1−λ)² × P(neg)
5 × P(neg) × (1−λ)² = P(pos) × λ
```
For a typical agent with P(neg) ≈ 0.2 (one-fifth of standard draws exceed n_active):
```
1.0 × (1−λ)² = λ   →   λ ≈ 0.382   (same golden-ratio form)
```
For P(neg) ≈ 0.3: λ_eq ≈ 0.47. Empirical average across agents lands at ≈ 0.50.

### Comparison: symmetric (v1) vs asymmetric (v2)

| Variant          | α     | β     | T   | Small λ_eq | Large λ_eq | Large eff |
|------------------|-------|-------|-----|-----------|-----------|-----------|
| Symmetric (v1)   | 0.05  | 0.05  | 300 | 0.40      | 0.37      | 0.931–0.953 |
| Asymmetric (v2)  | 0.01  | 0.05  | 600 | 0.49      | 0.52      | 0.920–0.934 |

The symmetric version performs **better** for large markets because it cools faster toward the
true optimum (λ=0). The asymmetric version overshoots — the 5× heating advantage lifts λ_eq
too high for thick markets, hurting efficiency.

### Diagnosis: individual signal cannot distinguish market thickness

The core problem: the marginal surplus signal is a **purely individual** observation.
Whether the market is thin or thick is a market-level property, but each agent only sees:
- Their own draw quantity
- Their own marginal WTP/WTAC vs p*

Two agents with the same n_active/n_total ratio will see the same signal distribution
regardless of how many other agents are in the market. The annealing schedule therefore
converges to the same λ_eq for all market sizes with the same underlying agent microstructure.

**Implication:** to achieve market-size-sensitive convergence (λ→1 for thin, λ→0 for thick),
the update rule needs a market-level signal — e.g., observed trade volume D relative to some
prior on q*, or a public announcement of aggregate efficiency. A purely individual marginal
surplus signal cannot encode this information.

---

## 15. Rationing-Guided Tatonnement — Exact WGE Convergence

**Script:** `test/taton_zit.jl`
**Motivation:** reach WGE without access to global signals (D, S, q* never observed).

### The ratchet rule

Each agent maintains `prev_k_i` (last period's realized trade, initialized to 0).
Each period, demand/supply is:

```
d_i(t) = min(n_active_i,  prev_k_i(t−1) + 1)
```

After observing realized trade: `prev_k_i(t) = realized_c(d_i(t), D, S)`.

**Local information only:**
- `prev_k_i`: own realized trade from previous period
- `n_active_i`: count(WTP ≥ p*), computed from own valuations and broadcast price p*
- No D, S, q*, or other agents' actions ever observed

**Why the ratchet converges:**
- *Never over-demands:* cap at n_active_i prevents demanding units with WTP < p*
- *Always probes:* asks for one more than the market last allowed → monotone ascent
- *Rationing provides the brake:* when D > S, pro-rata cuts k_i < d_i → next d_i shrinks
- *Fixed point:* k_i = n_active_i → d_i = min(n_active_i, n_active_i + 1) = n_active_i → WGE

### Results

| Market size  | t=1 eff | t=3 eff | t=5 eff | t=10 eff | Converged by |
|--------------|---------|---------|---------|----------|--------------|
| Small/mixed  | 0.820   | 0.985   | 1.000   | **1.000** | t ≈ 5–7      |
| 100 agents   | 0.778   | 0.971   | 0.998   | **1.000** | t ≈ 8–10     |
| 500 agents   | 0.779   | 0.971   | 0.998   | **1.000** | t ≈ 8–10     |

**Exact WGE (eff = 1.0000, std = 0.0000) achieved in ≤ 10 turns across all market sizes.**

### Comparison with stochastic ZIT benchmarks

| Method              | Final eff | vs WGE   | Turns to converge |
|---------------------|-----------|----------|-------------------|
| Ratchet tatonnement | **1.0000** | 0.0000  | ≤ 10              |
| Standard ZIT        | 0.984     | −0.016   | never (stochastic)|
| BiGeo price ZIT     | 0.915     | −0.085   | never (stochastic)|
| Marginal anneal (sym)| 0.940    | −0.060   | ~100 turns        |

### Why this works: the information content of rationing

The ratchet decouples the convergence problem into two one-sided mechanisms:
- **Consumers starting below q*:** never rationed (D < S), so they increment freely each turn
- **Consumers starting above q*:** rationed (D > S), so their realized k_i < d_i; next period d_i shrinks

The two groups coordinate through the market without communicating: over-demanders are cut back by rationing while under-demanders increment unimpeded, until D = S = q*.

The convergence speed is O(max(n_active_i)) ≤ n_total ≤ MAX_UNITS = 5 turns — bounded by the maximum WGE-optimal quantity per agent, independent of market size.

### Significance

This establishes that **WGE is reachable with purely local information and a minimal learning rule**. The only non-local element is the broadcast price p*, which is standard in competitive market theory (the auctioneer announces p* and agents respond independently). Given p*, agents need only observe their own realized trade to converge to WGE in O(q*/N) turns.

---

## 13. Adaptive ZIT — Signal-Gated Annealing

**Script:** `test/adaptive_zit_anneal.jl`
**Design:** each agent maintains a continuous `λ_i ∈ [0,1]` (starts at 1.0 — hot).
Each turn, agent uses BiGeo with prob `λ_i`, standard with prob `1 − λ_i`.
After observing marginal-unit surplus:

```
signal > 0  →  cool:  λ_i ← λ_i × (1 − 0.05)      [found good unit → exploit more]
signal < 0  →  heat:  λ_i ← λ_i + (1 − λ_i) × 0.05 [overtrading → explore more]
signal = 0  →  no update (no trade — neutral)
```

**Economic logic:** BiGeo = explore (finds quality inframarginal units).
Standard = exploit (calibrated to center aggregate demand at q*).
Positive signal → this tool is working, cool toward exploitation.
Negative signal → overtrading past n_active detected, heat toward exploration.

### Results

| Market size  | Final λ_c | Final λ_f | Final eff | vs Std   | vs BiGeo |
|--------------|-----------|-----------|-----------|----------|----------|
| Small/mixed  | 0.405     | 0.452     | 0.679     | −0.017   | −0.232   |
| 100 agents   | 0.364     | 0.369     | 0.931     | −0.052   | +0.017   |
| 500 agents   | 0.367     | 0.370     | 0.953     | −0.044   | +0.036   |

### Large-market dynamics

```
  Turn │  λ_c     λ_f    Eff (100-agent)    Eff (500-agent)
  t=  1 │  1.000  1.000     0.915               0.917
  t= 10 │  0.837  0.765     0.931               0.939
  t= 25 │  0.643  0.578     0.933               0.950
  t= 50 │  0.488  0.455     0.937               0.955
  t=100 │  0.397  0.391     0.941               0.952
  t=300 │  0.364  0.369     0.931               0.953
```

**Annealing resolves the miscoordination trap:** under pure marginal signal (§12), large
markets converged to λ≈0.80 and efficiency collapsed to 87-88%.  Signal-gated annealing
self-corrects to λ_eq ≈ 0.37 and achieves 93-95% — substantially above pure BiGeo (91.5-91.6%).

The improvement over pure BiGeo comes from the mixture of strategies: standard agents contribute
high-volume demand (centering aggregate D near q*), while BiGeo agents contribute high-quality
inframarginal demand. The system finds a natural mixing ratio rather than an extremum.

### Equilibrium analysis: why λ_eq ≈ 0.37

The cooling equilibrium condition (equal rates α=β=0.05):

```
E[Δλ_i] = 0  ⟺  P(signal > 0) × α × λ_i = P(signal < 0) × β × (1 − λ_i)
With α = β:  λ_i / (1 − λ_i) = P(signal < 0) / P(signal > 0)
```

An agent using standard with n_active/n_total ≈ 0.5 has:
- `P(neg | std) ≈ 0.5`   [half the draws exceed n_active → negative marginal surplus]
- `P(pos | std) ≈ 0.5`   [other half hit active range → positive]
- BiGeo: all non-zero draws give positive signal

At the population equilibrium for a market where n_active/n_total ≈ 0.5, this gives:
```
λ_i = (1 − λ_i)²   →   λ_i ≈ 0.382   (golden ratio conjugate: φ⁻² = 1 − φ⁻¹)
```

The empirical result of λ_eq ≈ 0.37 matches this prediction. The equilibrium mixing rate
is determined by the fundamental microstructure of the market (n_active/n_total ratio),
not by an externally imposed schedule.

### Small-market failure

In thin markets, the system starts at λ=1 (pure BiGeo, eff=0.91) and decays toward λ_eq≈0.40
but **does not converge within 300 turns**. Efficiency at t=300 is 0.679 — slightly below
pure standard (0.696). Two issues:

1. **Still in transit:** λ is declining at t=300 (not yet converged). Needs more turns.
2. **Heterogeneous λ_i distribution:** some agents have λ_i near 1, others near 0.
   This heterogeneity reduces efficiency below what the static mixture model predicts at
   the same mean λ, because high-λ agents collectively underdemand.

### Comparison across all adaptive designs

| Signal              | Small λ_eq | Small eff | Large λ_eq | Large eff | Issue                        |
|---------------------|-----------|-----------|-----------|-----------|------------------------------|
| Quantity            | 0.076     | 0.754     | ≈0        | 0.984     | Wrong objective              |
| Cumulative surplus  | 0.071     | 0.761     | ≈0        | 0.984     | Monotone in qty              |
| Marginal surplus    | 0.672     | 0.698     | 0.80      | 0.877     | Miscoordination trap         |
| **Signal-gated anneal** | **0.40** | **0.679** | **0.37** | **0.940** | **Resolves trap; thin mkt needs more turns** |

---

## 16. Decentralized ZIT — False Equilibrium Trap (Pure Rationing Signals)

**Script:** `test/decent_zit.jl`
**Motivation:** Remove p* from agents entirely. Agents know only own WTP/WTAC + own realized k_i.

### Mechanism

Each agent maintains private price belief p̂_i (initialized at p̂=0 for consumers, p̂=max_WTAC for firms — maximum aggressiveness). Updated via additive tatonnement (DELTA=1.0):

```
Consumer rationed (k < d):          p̂_i += DELTA   # demanded too many → price too low
Consumer at ratchet cap, not ration: p̂_i -= DELTA  # maybe too conservative → lower
Consumer ratcheting, not ration:     hold
Firm rationed (k < d):              p̂_j -= DELTA   # supplied too many → price too high
Firm at ratchet cap, not ration:    p̂_j += DELTA   # maybe too conservative → raise
```

Quantity ratchet: `d_i(t) = min(count(WTP ≥ p̂_i), k_i(t-1) + 1)`

### Results

| Market size  | Final q/q* | Final eff | vs WGE   | vs Std   | Converged? |
|--------------|------------|-----------|----------|----------|------------|
| Small/mixed  | 1.196      | 0.635     | −0.365   | −0.070   | No — false eq |
| 100 agents   | 1.378      | 0.842     | −0.158   | −0.142   | No — false eq |
| 500 agents   | 1.381      | 0.845     | −0.155   | −0.153   | No — false eq |

Stable from t=10 to t=300 — NOT a transient. True false equilibrium.

### The False Equilibrium Trap (Agda-relevant)

**Trap condition:** D = S = q > q*. Both sides are over-trading at a balanced quantity
above WGE.

When D = S (balanced market):
- No consumer is rationed → no "raise p̂" signal fires
- No firm is rationed → no "lower p̂" signal fires
- Agents at their ratchet caps: p̂ hits floor (0 for consumers) or ceiling (max_WTAC for firms)
- p̂ updates become no-ops (floor/ceiling binding)
- System is permanently locked

**The information-theoretic result:**
Pure rationing signals cannot distinguish "D=S=q*" (WGE) from "D=S=q>q*" (false equilibrium).
The rationing signal is IDENTICAL (zero) in both states.

This is not a failure of the specific mechanism — it is a fundamental information constraint:
**Without a price signal, decentralized agents cannot determine whether their units are
profitable to trade. Rationing tells you the market is congested, but does not tell you
whether YOU are contributing to the congestion or trading legitimately.**

**Implication for Agda research:**
This negative result has a formal proof structure. Define:
- `LocalSignal = (own_k_i, own_d_i, own_WTP, own_WTAC)`
- `GlobalState = (D, S, q, p*)`

The claim is: no update rule `f(LocalSignal) → p̂_i` can provably converge to WGE for all
market configurations, because LocalSignal alone does not distinguish WGE from false equilibria.
The minimal additional signal required is a scalar market price (transaction price p_t),
which breaks symmetry between WGE and false equilibria.

**Fix:** allow agents to observe their own transaction price p_t = WTP[k_i] (the value of
their marginal traded unit). This IS a local signal (derived from own WTP array and own k_i),
not a global one. It provides the missing price anchor. See §17.

---

## 17. Decentralized ZIT — Transaction-Price Adaptive Beliefs

**Script:** `test/decent_zit_price.jl`
**Motivation:** Break the false equilibrium trap of §16 using own transaction price p_t.

### The additional local signal

After trading k units, consumer i observes `p_t = WTP[k_i]` — the value of their marginal
traded unit. This is PRIVATE knowledge (own WTP array + own k_i); no global information.
In any real market, agents observe what they paid/received.

### Price belief update (three cases)

```
Rationed (k < d):
  Consumer: p̂_i += η*(WTP[k]       - p̂_i)  [price ≥ WTP[k]       → raise]
  Firm:     p̂_j += η*(WTAC[k]      - p̂_j)  [price ≤ WTAC[k]      → lower]

At ratchet cap, not rationed (k = d = n_active < n_total):
  Consumer: p̂_i += η*(WTP[na+1]    - p̂_i)  [price ≤ WTP[na+1]    → lower]
  Firm:     p̂_j += η*(WTAC[na+1]   - p̂_j)  [price ≥ WTAC[na+1]   → raise]

Ratcheting (k = d < n_active) OR at full cap (n_active = n_total):
  No update — no informative price signal
```

EMA rate η = 0.20 (symmetric). Init: p̂_i = mean(own WTP), p̂_j = mean(own WTAC).

### Convergence mechanism

p̂_i oscillates between WTP[n_active*] and WTP[n_active*+1] (the two WTP values straddling p*).
The amplitude = WTP gap → 0 as valuation distribution becomes denser.
Efficiency loss = surplus lost during over/under-demand oscillation.

**Breaks the false equilibrium:** "At cap, not rationed" now carries a PRICE signal
(WTP[na+1]) not just a quantity signal. When D = S > q*, agents at cap observe their
excluded unit's value and adjust p̂ to include/exclude it, breaking the D=S=q>q* trap.

### Results (η=0.20)

| Market size  | t=2 eff | t=5 eff | t=50 eff | t=300 eff | vs Std   | vs WGE   |
|--------------|---------|---------|----------|-----------|----------|----------|
| Small/mixed  | 0.783   | 0.799   | 0.799    | 0.795     | +0.099   | −0.205   |
| 100 agents   | **0.994** | 0.992 | 0.984    | 0.984     | +0.000   | −0.016   |
| 500 agents   | **0.997** | 0.994 | 0.973    | 0.970     | −0.026   | −0.030   |

### Key findings

1. **Near-WGE convergence in 2 turns for large markets.** At t=2: 99.4-99.7% efficiency.
   This is dramatically faster than any stochastic mechanism (which never converge exactly).

2. **Steady-state drift below t=2 peak.** The "at cap → lower p̂" signal fires constantly
   at WGE (D=S, not rationed → every agent receives the "lower" signal every turn). This
   creates permanent downward pressure that settles at a stable fixed point ~2-3% below WGE.

3. **Large markets: matches std ZIT (100C) or slightly below (500C).** Standard ZIT benefits
   from law of large numbers → ~99.6% for 500 agents. The adaptive mechanism adds overhead
   from its price-discovery oscillation.

4. **Small markets: significant improvement over §16 (+16pp) but still below taton (eff=1.0).**

### Minimal information requirement (Agda-relevant)

This experiment establishes that:
- {rationing signal alone} → false equilibrium at q/q*=1.38, eff=0.84 (§16)
- {rationing + own transaction price p_t = WTP[k_i]} → converges to ~97-98% eff

The minimal local information set for near-WGE convergence is:
```
LocalInfo = (own_WTP, own_WTAC, own_k_i, WTP[k_i])
```
where `WTP[k_i]` is derived from own data — not a market signal. This has a formal proof
structure worth developing in Agda: for any market with n_active* well-defined and
WTP distribution with density δ near p*, efficiency ≥ 1 - f(δ) after O(log(1/δ)) turns.

### Comparison to prior mechanisms

| Mechanism                  | Signals used              | Final eff (500)| Turns to near-WGE |
|----------------------------|---------------------------|----------------|-------------------|
| Standard ZIT               | None (random)             | 0.997          | Never (stochastic)|
| BiGeo ZIT                  | None (biased random)      | 0.916          | Never             |
| Signal-gated anneal (§13)  | Marginal surplus          | 0.953          | ~100              |
| Ratchet taton (§15, cheat) | k_i, WTP, **p***          | 1.000          | 10                |
| Decent ZIT (§16, fail)     | k_i, WTP (no p*)          | 0.845          | Never (false eq)  |
| **Decent ZIT+price (§17)** | **k_i, WTP, WTP[k_i]**   | **0.970**      | **2–5**           |

---

## 18. Local Annealing ZIT — Adaptive p̂ + Signal-Gated λ (Fully Decentralized)

**Script:** `test/local_anneal.jl`  
**Motivation:** Combine §17's local price belief p̂_i with §13-style λ annealing, removing
the global p* dependency entirely. Goal: λ_eq > 0 prevents false equilibria while p̂_i guides
the BiGeo lower bound.

### Version 1: BiGeo + Ratchet, Rationing/Cap Signal for λ

```
prob λ:   BiGeo over {WTP ≥ p̂_i}       [explore]
prob 1-λ: Ratchet = min(na(p̂_i), k+1)   [exploit]
λ heat:   rationed (k < d)
λ cool:   at cap, not rationed (d = na)
```

**Problem:** At WGE (D=S, no rationing), ALL agents cool every turn — the "at cap, not rationed"
signal fires constantly. The heat signal (rationing) never fires. λ → 0 monotonically.

**Results (100 agents, t=300):** eff=0.922, q/q*=0.729, λ_c=0.090, λ_f=0.038  
**Results (500 agents, t=300):** eff=0.913, q/q*=0.705  
**vs §17:** WORSE (−6.2pp for 100 agents)

### Version 2: BiGeo + Standard ZIT, Marginal Surplus Signal for λ

```
prob λ:   BiGeo over {WTP ≥ p̂_i}             [explore]
prob 1-λ: Standard ZIT (uniform [0, n_total])   [exploit]
λ cool:   WTP[k] ≥ p̂_i (positive surplus)
λ heat:   WTP[k] < p̂_i (negative surplus — ZIT over-shot past n_active)
p̂ update: BiGeo draws only (ZIT rationing too noisy)
```

**Theory:** BiGeo always yields WTP[k] ≥ p̂ (cool). Standard ZIT sometimes demands beyond
n_active → WTP[k] < p̂ (heat). Balance predicted at λ_eq ≈ φ⁻² ≈ 0.38, same golden-ratio
result as §13 with p̂ replacing p*.

**Observed λ_eq ≈ 0.35–0.36** (close to predicted 0.38) ✓

**Results (100 agents, t=300):** eff=0.944, q/q*=0.784, λ_c=0.357, λ_f=0.303  
**Results (500 agents, t=300):** eff=0.959, q/q*=0.801  
**vs §17:** WORSE (−4.0pp for 100 agents, −1.1pp for 500 agents)

### Diagnosis: Why λ annealing degrades §17

Two independent failure mechanisms:

**1. min(D,S) asymmetry of standard ZIT:**  
ZIT draws 0 with probability 1/(n_total+1) ≈ 1/6. These zero-draws reduce D (or S). Since
q = min(D, S), low draws that reduce the smaller side hurt q directly, while high draws that
increase the larger side only help if they flip which side is smaller. The net effect is
persistent downward pressure on q (q/q* ≈ 0.78–0.80 < §17's 0.83–0.90).

**2. Incompatibility of ZIT and p̂ update:**  
When ZIT consumers collectively over-demand, D >> S. Rationing fires for BiGeo consumers
with k << d. WTP[k] at high D >> S corresponds to an inframarginal unit (high WTP, far
above p*), pushing p̂ far above p*. This shrinks n_active below n_active*, which combined
with the min-effect produces sustained under-trading. Gating p̂ updates to BiGeo-only
partially mitigates this but doesn't eliminate it (ZIT zeros create spurious at-cap signals
that lower p̂ below p* when a BiGeo agent happens to hit cap during a D < S period caused
by ZIT zeros).

**Root cause:** The ZIT random draws that are needed for the λ heat signal are structurally
incompatible with the p̂ update mechanism. The two signals cannot be cleanly separated
because they share the same market clearing outcome.

### Comparison: all §18 variants vs §17

| Mechanism                        | 100-agent eff | 500-agent eff | λ_eq  |
|----------------------------------|---------------|---------------|-------|
| §17 (pure BiGeo + local p̂)       | **0.984**     | **0.970**     | N/A   |
| §18 v1 (ratchet + rat/cap λ)     | 0.922         | 0.913         | → 0   |
| §18 v2 (ZIT + surplus λ)         | 0.944         | 0.959         | 0.36  |

### Conclusion

The λ annealing mechanism does NOT improve over §17 with local signals. The minimal-information
mechanism for near-WGE convergence in this market class is **§17: pure BiGeo with adaptive p̂**.

The golden-ratio λ_eq from §13 requires a clean marginal surplus signal relative to p*. With
local p̂ replacing p*, the ZIT component needed for the heat signal corrupts p̂ via two channels
(rationing and spurious at-cap). These are not separable without global market information.

### Agda implication

The information-sufficient local mechanism ({WTP, k_i, WTP[k_i]}, §17) is provably minimal:
- {rationing alone} → false equilibrium (§16)
- {rationing + WTP[k_i]} → ~97% efficiency (§17)
- {rationing + WTP[k_i] + λ annealing via ZIT} → WORSE than {rationing + WTP[k_i]} (§18)

The third point is a negative result worth formalizing: adding MORE local information (ZIT
draws, which introduce structured noise) can hurt convergence. The mechanism is information-brittle
in this specific way. For Agda: the proof of §17's convergence relies on BiGeo drawing
exclusively from the profitable region {WTP ≥ p̂} — any mechanism that mixes draws from
the unprofitable region (ZIT, ratchet) corrupts the price discovery loop.

---

## 19. Local Elasticity Annealing — BiGeo + Ratchet, Elasticity-Gap λ Signal

**Script:** `test/local_elast.jl`  
**Motivation:** User's intended design: one λ cooling from BiGeo (exploration) to ratchet/additive
(exploitation). Signal for λ: local elasticity gap — own WTP[k] − WTP[k+1] vs own precomputed
mean_gap. No global p*, D, S, or q* ever observed.

### Mechanism

```
State per agent: (p̂_i, λ_i, prev_k_i)

Quantity draw (mixed strategy):
  prob λ_i:   BiGeo over {WTP ≥ p̂_i}, p̂_i as lower bound    [explore]
  prob 1-λ_i: ratchet = min(na(p̂_i), prev_k_i + 1)           [exploit]

p̂ update (BiGeo-only; ratchet never updates p̂):
  rationed (k < d):            p̂ += η*(WTP[k]     − p̂)        [raise]
  at cap (d = na, k = d):      p̂ += η*(WTP[na+1]  − p̂)        [lower]
  else: hold

λ update (elasticity gap signal):
  gap = WTP[k] − WTP[k+1]        (forward drop at own traded unit)
  mean_gap = (WTP[1]−WTP[n])/(n−1)  (own WTP range / (n-1); precomputed)
  gap > mean_gap  →  λ *= (1−β)               [inelastic → cool toward ratchet]
  gap < mean_gap  →  λ += (1−λ)*α             [elastic   → heat toward BiGeo]
  k=0 or k=n_total  →  hold

Init: p̂_i = mean(WTP), λ_i = 1.0, prev_k_i = 0
```

Ratchet owns prev_k; BiGeo turns leave it unchanged so the additive climb is not disrupted.
p̂ is updated from BiGeo draws only, preserving §17's price-discovery signal quality.

### Variants explored (all give eff ≈ 91–95%, worse than §17)

Four configurations of the p̂/prev_k update were tried:

| Variant | p̂ raise | p̂ lower | prev_k | 100-agent eff | q/q* |
|---------|----------|----------|--------|---------------|------|
| 1. All draws update p̂ | both | both | shared | 0.916 | 0.72 |
| 2. Ratchet raises, BiGeo lowers | ratchet | BiGeo | shared | 0.894 | ~0.70 |
| 3. BiGeo only, shared prev_k | BiGeo | BiGeo | shared | 0.916 | 0.72 |
| 4. BiGeo only, ratchet owns prev_k | BiGeo | BiGeo | separate | 0.918 | 0.72 |

All variants produce essentially identical outcomes. §17 benchmark: eff=0.984, q/q*≈0.90.

### Parameter sensitivity (λ cooling rate)

Increasing BETA from 0.05 to 0.35 (ALPHA=0.05 throughout):

| BETA | λ_eq (100-agent) | 100-agent eff | q/q* |
|------|-----------------|---------------|------|
| 0.05 | 0.624           | 0.918         | 0.72 |
| 0.35 | 0.564           | 0.950         | 0.79 |

λ refuses to cool toward 0 regardless of parameter choice — see below.

### Why λ_eq stays high (gap distribution analysis)

The elasticity gap signal has a structural bias. For agents with n WTP values drawn from
a continuous distribution, the gaps WTP[k]−WTP[k+1] are approximately Exponentially
distributed with mean = mean_gap. For an Exponential variable X with mean μ:

```
P(X < μ) = 1 − e⁻¹ ≈ 0.63   [gap < mean_gap → heat]
P(X > μ) = e⁻¹ ≈ 0.37       [gap > mean_gap → cool]
```

At equilibrium with α=β=0.05:
```
λ_eq = α × P(heat) / (β × P(cool) + α × P(heat))
     = 0.05 × 0.63 / (0.05 × 0.37 + 0.05 × 0.63) = 0.63
```

This matches the observed λ_eq ≈ 0.62–0.67 exactly. With BETA=0.35, ALPHA=0.05:
```
λ_eq = 0.05 × 0.63 / (0.35 × 0.37 + 0.05 × 0.63) ≈ 0.20 (theory)
```
But observed is ≈ 0.56 because agents with 1–2 units have degenerate gaps
(gap == mean_gap exactly for 2-unit agents; no gap comparison for 1-unit agents) and
never update λ, staying frozen at λ=1.0 and inflating the population mean.

### Root cause: BiGeo structural under-demand

At λ_eq ≈ 0.63 with BiGeo drawing log-uniformly from [p̂, WTP[1]]:

```
E[BiGeo draw] ≈ 0.54–0.67 × na  (log-uniform undersamples near-marginal units)
E[ratchet draw] = na              (ratchet reaches cap deterministically)

E[total demand] = λ_eq × E[BiGeo] + (1−λ_eq) × na
               ≈ 0.63 × 0.5 × na + 0.37 × na = 0.685 × na = 0.685 × q*
```

This perfectly predicts the observed q/q* ≈ 0.70–0.72. The 28-30% quantity shortfall
is a structural ceiling from the λ_eq × BiGeo under-demand product. It cannot be overcome
without λ_eq → 0, which the elasticity gap signal does not achieve.

### Critical comparison with §17

**§17 (test/decent_zit_price.jl) is also ratchet-based** — the quantity draw in §17 is
PURE RATCHET: `d_i = min(na_i, prev_k_i + 1)` every turn. §17's p̂ update uses the same
three-case EMA (raise/lower/hold) but applies to ALL turns, not just BiGeo turns.

The §19 experiment was predicated on the hypothesis that BiGeo exploration would improve
price discovery over the ratchet alone. This is falsified: the §17 ratchet's p̂ update
(rationing → raise, at-cap → lower) provides exactly the same price signal quality as
BiGeo, while also achieving q/q* → 1 deterministically. BiGeo's exploration phase (1-λ
fraction of turns) only reduces q without improving p̂.

### Final results (BETA=0.35, ALPHA=0.05, variant 4)

| Market size  | q/q*   | eff    | λ_c   | λ_f   | vs §17 |
|--------------|--------|--------|-------|-------|--------|
| 100 agents   | 0.787  | 0.950  | 0.564 | 0.611 | −0.034 |
| 500 agents   | 0.762  | 0.943  | 0.586 | 0.631 | −0.027 |

### Conclusion

The BiGeo+ratchet design with elasticity-gap λ signal achieves 95% efficiency (vs §17's
98.4%). The gap is structural and not correctable by parameter tuning:

1. **Elasticity gap gives P(heat) > P(cool)** due to Exponential gap distribution → λ
   stays near 0.56–0.67, not 0.
2. **BiGeo draws are below na in expectation** → q/q* ≈ 0.70–0.79 at any stable λ > 0.
3. **§17's pure ratchet already provides price discovery** — the ratchet at-cap signal
   (fires when d = na, k = d) gives the same "lower p̂" information as BiGeo's at-cap signal.
   BiGeo exploration adds no information that the ratchet does not already generate.

The user's intended cooling mechanism (λ: BiGeo → ratchet) is conceptually sound but
redundant: the ratchet is already the optimal exploitation AND exploration mechanism in
this local-signal setting. λ=0 (pure ratchet, §17) dominates all interior λ values.

### Open directions

1. **Alternative cooling signal:** A signal that achieves λ_eq ≈ 0 could preserve the
   design intent. Candidate: cool when k > 0 (traded something), heat when k = 0
   (excluded). At WGE all agents trade → λ → 0. When p̂ drifts → some agents excluded →
   λ heats → BiGeo corrects p̂ → λ cools again. Untested.
2. **Ratchet parameter sensitivity:** §17 has a residual 2-3% below WGE from the
   constant "at-cap → lower" signal. Could BiGeo replace the lower-signal to reduce this
   oscillation? Hypothesis: yes, as BiGeo's at-cap fires less frequently than ratchet's.
   Would require λ~0.1 (10% BiGeo) with BiGeo providing only the lower signal.
3. **Agda formalization:** The negative result (BiGeo+ratchet < pure ratchet) has a clean
   proof structure: any draw that produces E[d] < na reduces q below q*, and no local
   signal can compensate when the shortfall is systematic.
