# BiGeometric ZIT — Price-Space Formulation.
#
# In standard ZIT, a consumer draws quantity uniformly from [0, Q_max].
# Here we ask: what is the BiGeometric analog in PRICE space?
#
# Each consumer draws a "virtual threshold price" p_v log-uniformly from
# [p*, max_WTP_i], then demands all units with WTP ≥ p_v.
#
# In log-price space this is just Uniform[log(p*), log(max_wtp)].
# The resulting probability of buying unit j (with wtp_j ≥ p*) is:
#
#   P(buy unit j) = log(wtp_j / p*) / log(wtp_1 / p*)
#
# where wtp_1 = max_wtp_i is the consumer's highest valuation.
# Unit 1 (highest WTP) is ALWAYS bought (P = 1).
# Units closer to p* in multiplicative terms are bought less often.
#
# Key structural property:
#   BiGeo demand ≤ WGE demand always  →  efficiency ∈ [0, 1] always.
#   No negative efficiency is possible.
#
# Firms draw p_v ~ LogUniform[min_WTAC_i, p*] and supply all units with
# WTAC ≤ p_v, giving P(supply unit j) = log(p* / wtac_j) / log(p* / wtac_1).
#
# Run with:
#   julia --project=. test/bgi_price_efficiency.jl

using Random
using Statistics
using Printf
using DiscreteMarket

const N_MARKETS = 1_000
const T_SIM     = 1_000
const SEED_BASE = 20_000    # same seeds as zi_efficiency_analysis.jl

# ── Surplus helper ────────────────────────────────────────────────────────────

function surplus_f64(m::GoodMarket, q::Int) :: Float64
    q <= 0 && return 0.0
    all_wtp  = sort([Float64(v) for d in m.consumers for v in d.wtp], rev=true)
    all_wtac = sort([Float64(c) for f in m.firms     for c in f.wtac])
    n = min(q, length(all_wtp), length(all_wtac))
    n == 0 && return 0.0
    sum(all_wtp[1:n]) - sum(all_wtac[1:n])
end

# ── BiGeometric price-space demand for one consumer ───────────────────────────
# Active units: wtp ≥ p*.  Draw p_v ~ LogUniform[p*, max_active_wtp].
# Return count of active units with wtp ≥ p_v.
# Edge case: if all active units share the same wtp (= max), p_v = p* always
# and demand = all active units deterministically.

function bgi_price_consumer_demand(d::ConsumerDemand, p_star::Price,
                                    rng::AbstractRNG) :: Int
    # Collect active units (WTP ≥ p*)
    active = [Float64(v) for v in d.wtp if v >= p_star]
    isempty(active) && return 0

    max_wtp = active[1]             # wtp is nonincreasing, so first = max
    log_p   = log(Float64(p_star))
    log_max = log(max_wtp)

    if log_max ≈ log_p              # all active units at exactly p* — no spread
        return length(active)
    end

    u     = rand(rng)
    log_pv = log_p + u * (log_max - log_p)
    pv    = exp(log_pv)

    count(w -> w >= pv, active)
end

# ── BiGeometric price-space supply for one firm ───────────────────────────────
# Active units: wtac ≤ p*.  Draw p_v ~ LogUniform[min_active_wtac, p*].
# Return count of active units with wtac ≤ p_v.

function bgi_price_firm_supply(f::FirmSupply, p_star::Price,
                                rng::AbstractRNG) :: Int
    active = [Float64(c) for c in f.wtac if c <= p_star]
    isempty(active) && return 0

    min_wtac = active[1]            # wtac is nondecreasing, so first = min
    log_p    = log(Float64(p_star))
    log_min  = log(min_wtac)

    if log_p ≈ log_min              # all active units at exactly p*
        return length(active)
    end

    u      = rand(rng)
    log_pv  = log_min + u * (log_p - log_min)
    pv     = exp(log_pv)

    count(c -> c <= pv, active)
end

# ── Simulate T periods under BiGeo price-space ZIT ───────────────────────────

function bgi_price_simulate(m::GoodMarket, p_star::Price, T::Int;
                             seed::Int=0) :: Vector{Int}
    rng    = MersenneTwister(seed)
    traded = Vector{Int}(undef, T)
    for t in 1:T
        demand = sum(bgi_price_consumer_demand(d, p_star, rng) for d in m.consumers;
                     init=0)
        supply = sum(bgi_price_firm_supply(f, p_star, rng) for f in m.firms;
                     init=0)
        traded[t] = min(demand, supply)
    end
    traded
end

# ── Simulate T periods under standard ZIT (for head-to-head on same seed) ────

function std_simulate(m::GoodMarket, p_star::Price, T::Int; seed::Int=0) :: Vector{Int}
    zi_consumers = [ZIConsumer(1, p_star * length(d.wtp)) for d in m.consumers]
    zi_firms     = [ZIFirm(1, length(f.wtac))             for f in m.firms]
    zm  = ZIMarket(zi_consumers, zi_firms, 1)
    res = zi_simulate(zm, [p_star], T; seed=seed)
    res.traded[:, 1]
end

# ── Theoretical per-unit purchase probability ─────────────────────────────────
# P(buy unit j) = log(wtp_j / p*) / log(wtp_1 / p*)
# Expected demand = Σ_j P(buy unit j) over active units.

function bgi_price_expected_demand(d::ConsumerDemand, p_star::Price) :: Float64
    active = [Float64(v) for v in d.wtp if v >= p_star]
    isempty(active) && return 0.0
    max_wtp = active[1]
    log_range = log(max_wtp) - log(Float64(p_star))
    log_range ≈ 0.0 && return Float64(length(active))
    sum(log(w) - log(Float64(p_star)) for w in active) / log_range
end

function bgi_price_expected_supply(f::FirmSupply, p_star::Price) :: Float64
    active = [Float64(c) for c in f.wtac if c <= p_star]
    isempty(active) && return 0.0
    min_wtac = active[1]
    log_range = log(Float64(p_star)) - log(min_wtac)
    log_range ≈ 0.0 && return Float64(length(active))
    sum(log(Float64(p_star)) - log(c) for c in active) / log_range
end

# ── Record ────────────────────────────────────────────────────────────────────

struct PriceRecord
    q_wge       :: Int
    eff_std     :: Float64
    eff_bgi_p   :: Float64
    mean_q_std  :: Float64
    mean_q_bgip :: Float64
    expected_d  :: Float64          # theoretical E[demand] for BiGeo price
    expected_s  :: Float64          # theoretical E[supply] for BiGeo price
end

# ── Main ──────────────────────────────────────────────────────────────────────

function run_price_comparison()

records    = PriceRecord[]
n_no_equil = 0
n_zero_s   = 0

for seed in 1:N_MARKETS
    rng   = MersenneTwister(seed + SEED_BASE)
    n_c   = rand(rng, 1:8)
    n_f   = rand(rng, 1:8)
    max_u = rand(rng, 1:6)

    m, r = generate_good_market(rng;
        good=1, n_consumers=n_c, n_firms=n_f, max_units=max_u, Q=100)

    r.cleared || (n_no_equil += 1; continue)
    p_star = r.price

    q_star = aggregate_demand(m, p_star)
    s_wge  = surplus_f64(m, q_star)
    s_wge > 0 || (n_zero_s += 1; continue)

    # Standard ZIT
    traded_std = std_simulate(m, p_star, T_SIM; seed=seed)
    eff_std    = mean(surplus_f64(m, t) for t in traded_std) / s_wge

    # BiGeo price-space ZIT
    traded_bgip = bgi_price_simulate(m, p_star, T_SIM; seed=seed)
    eff_bgip    = mean(surplus_f64(m, t) for t in traded_bgip) / s_wge

    # Theoretical expected quantities
    exp_d = sum(bgi_price_expected_demand(d, p_star) for d in m.consumers)
    exp_s = sum(bgi_price_expected_supply(f, p_star) for f in m.firms)

    (isnan(eff_std) || isnan(eff_bgip)) && continue

    push!(records, PriceRecord(
        q_star, eff_std, eff_bgip,
        mean(Float64.(traded_std)),
        mean(Float64.(traded_bgip)),
        exp_d, exp_s
    ))
end

n = length(records)
effs_std  = [r.eff_std   for r in records]
effs_bgip = [r.eff_bgi_p for r in records]

println("=" ^ 64)
println("BiGeometric ZIT — Price-Space Formulation")
println("  Consumer: p_v ~ LogUniform[p*, max_WTP],  buy if wtp ≥ p_v")
println("  Firm:     p_v ~ LogUniform[min_WTAC, p*], sell if wtac ≤ p_v")
println("=" ^ 64)
println()

# ── Theoretical demand illustration ──────────────────────────────────────────

println("── P(buy unit j) = log(wtp_j/p*) / log(wtp_max/p*) ─────────")
println("  Example: consumer with wtp = [4/1, 3/1, 2/1] at p* = 1/1:")
println("    Unit 1 (wtp=4): P = log(4)/log(4) = 1.000")
println("    Unit 2 (wtp=3): P = log(3)/log(4) = 0.792")
println("    Unit 3 (wtp=2): P = log(2)/log(4) = 0.500")
println("    E[demand] = 2.292  vs WGE demand = 3, Standard ZIT mean = 1.5")
println()

# ── Side-by-side efficiency summary ──────────────────────────────────────────

println("── Efficiency: Standard vs BiGeo Price-Space (n=$n markets) ──")
println()
println("                Standard ZIT    BiGeo Price     Δ")
@printf("  Mean:         %.4f          %.4f        %+.4f\n",
        mean(effs_std), mean(effs_bgip), mean(effs_bgip) - mean(effs_std))
@printf("  Median:       %.4f          %.4f        %+.4f\n",
        median(effs_std), median(effs_bgip), median(effs_bgip) - median(effs_std))
@printf("  Std:          %.4f          %.4f        %+.4f\n",
        std(effs_std), std(effs_bgip), std(effs_bgip) - std(effs_std))
@printf("  Min:          %.4f         %.4f       %+.4f\n",
        minimum(effs_std), minimum(effs_bgip), minimum(effs_bgip) - minimum(effs_std))
@printf("  Max:          %.4f          %.4f        %+.4f\n",
        maximum(effs_std), maximum(effs_bgip), maximum(effs_bgip) - maximum(effs_std))
println()
@printf("  Mean/Std (Sharpe):  Std = %.3f,  BiGeo-P = %.3f\n",
        mean(effs_std)/std(effs_std), mean(effs_bgip)/std(effs_bgip))
println()

println("── Bounds check (efficiency must be ∈ [0,1] for BiGeo price) ─")
n_neg = count(e -> e < -1e-9, effs_bgip)
n_over = count(e -> e > 1 + 1e-9, effs_bgip)
@printf("  Negative efficiency: %d / %d\n", n_neg, n)
@printf("  Efficiency > 1:      %d / %d\n", n_over, n)
println()

println("── Tail counts ───────────────────────────────────────────────")
println("  Threshold │  Standard    BiGeo-Price")
for thresh in [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    n_std = count(e -> e >= thresh, effs_std)
    n_bgip = count(e -> e >= thresh, effs_bgip)
    @printf("  eff≥%.2f  │  %3d/%d (%.1f%%)  %3d/%d (%.1f%%)\n",
            thresh,
            n_std, n, 100*n_std/n,
            n_bgip, n, 100*n_bgip/n)
end
println()

# ── Breakdown by q* ──────────────────────────────────────────────────────────

println("── By WGE quantity q* ────────────────────────────────────────")
println("  q*     │  n    Std mean   BiGeo-P mean   Δ")
for (rng_q, label) in [(1:2,"1–2"),(3:5,"3–5"),(6:10,"6–10"),(11:typemax(Int),"11+")]
    sub = [r for r in records if r.q_wge in rng_q]
    isempty(sub) && continue
    ms = mean(r.eff_std   for r in sub)
    mb = mean(r.eff_bgi_p for r in sub)
    @printf("  q*=%-5s │  %-4d  %.4f     %.4f        %+.4f\n",
            label, length(sub), ms, mb, mb - ms)
end
println()

# ── Head-to-head ─────────────────────────────────────────────────────────────

n_bgip_wins = count(r -> r.eff_bgi_p > r.eff_std, records)
n_std_wins  = count(r -> r.eff_bgi_p < r.eff_std, records)
n_tie       = count(r -> r.eff_bgi_p ≈ r.eff_std, records)
println("── Head-to-head (same market, same seed) ─────────────────────")
@printf("  BiGeo-P > Std:  %d / %d  (%.1f%%)\n", n_bgip_wins, n, 100*n_bgip_wins/n)
@printf("  Std > BiGeo-P:  %d / %d  (%.1f%%)\n", n_std_wins,  n, 100*n_std_wins/n)
@printf("  Tie:            %d / %d  (%.1f%%)\n", n_tie,       n, 100*n_tie/n)
println()

# ── Quantity comparison ───────────────────────────────────────────────────────

q_wge_all   = Float64.([r.q_wge       for r in records])
q_std_all   = [r.mean_q_std  for r in records]
q_bgip_all  = [r.mean_q_bgip for r in records]
exp_d_all   = [r.expected_d  for r in records]
exp_s_all   = [r.expected_s  for r in records]

println("── Quantity means vs WGE q* ──────────────────────────────────")
@printf("  WGE q*:            mean=%.2f\n", mean(q_wge_all))
@printf("  Standard ZIT q:    mean=%.2f  (%.1f%% of q*)\n",
        mean(q_std_all),  100*mean(q_std_all) /mean(q_wge_all))
@printf("  BiGeo-P traded q:  mean=%.2f  (%.1f%% of q*)\n",
        mean(q_bgip_all), 100*mean(q_bgip_all)/mean(q_wge_all))
@printf("  BiGeo-P E[demand]: mean=%.2f  (theoretical)\n", mean(exp_d_all))
@printf("  BiGeo-P E[supply]: mean=%.2f  (theoretical)\n", mean(exp_s_all))
println()

# ── WTP heterogeneity analysis ────────────────────────────────────────────────
# BiGeo price advantage comes from respecting WTP spread.
# Measure spread via coefficient of variation of WTP values.

println("── WTP heterogeneity and BiGeo-P efficiency ──────────────────")
println("  (How does WTP spread within a consumer affect BiGeo-P behavior?)")
println()

# Compute mean log-spread per market: log(max_wtp / p*) for active consumers
log_spreads = Float64[]
for seed in 1:N_MARKETS
    rng   = MersenneTwister(seed + SEED_BASE)
    n_c   = rand(rng, 1:8)
    n_f   = rand(rng, 1:8)
    max_u = rand(rng, 1:6)
    m, r  = generate_good_market(rng; good=1, n_consumers=n_c, n_firms=n_f,
                                  max_units=max_u, Q=100)
    r.cleared || continue
    p     = r.price
    surplus_f64(m, aggregate_demand(m, p)) > 0 || continue
    for d in m.consumers
        active = [Float64(v) for v in d.wtp if v >= p]
        (isempty(active) || length(active) <= 1) && continue
        push!(log_spreads, log(active[1]) - log(Float64(p)))
    end
end

if !isempty(log_spreads)
    println("  Log-spread log(max_wtp/p*) across all active consumers:")
    @printf("    Mean:   %.3f  (ratio = %.1fx)\n",
            mean(log_spreads), exp(mean(log_spreads)))
    @printf("    Median: %.3f  (ratio = %.1fx)\n",
            median(log_spreads), exp(median(log_spreads)))
    @printf("    P(spread ≈ 0):  %.1f%%  (homogeneous WTP → P(buy)=1 for all units)\n",
            100 * count(s -> s < 0.01, log_spreads) / length(log_spreads))
end
println()

# ── Three-way summary ─────────────────────────────────────────────────────────

println("=" ^ 64)
println("Three-way summary (same 1000 markets)")
println("=" ^ 64)
println()
println("  Method          Mean    Median   Std    Min      Sharpe")
println("  ─────────────────────────────────────────────────────────")
# Standard ZIT (from our run)
@printf("  Standard ZIT    %.4f  %.4f  %.4f  %7.4f  %.3f\n",
        mean(effs_std), median(effs_std), std(effs_std),
        minimum(effs_std), mean(effs_std)/std(effs_std))
# BiGeo price-space
@printf("  BiGeo Price     %.4f  %.4f  %.4f  %7.4f  %.3f\n",
        mean(effs_bgip), median(effs_bgip), std(effs_bgip),
        minimum(effs_bgip), mean(effs_bgip)/std(effs_bgip))
println()
println("  (BiGeo Quantity results from bgi_efficiency.jl for reference:)")
println("  BiGeo Quantity  0.5446  0.5676  0.2416  -5.1830  2.254")
println()

end # run_price_comparison

run_price_comparison()
