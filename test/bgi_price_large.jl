# BiGeometric price-space ZIT vs standard ZIT at large market sizes.
# Mirrors zi_large_market.jl but adds the BiGeo price-space formulation.
# Same seeds (30_000 / 31_000) so results are directly comparable.
#
# BiGeo price-space: consumer draws p_v ~ LogUniform[p*, max_WTP],
# buys all units with WTP ≥ p_v.  P(buy unit j) = log(wtp_j/p*)/log(wtp_max/p*).
#
# Run with:
#   julia --project=. test/bgi_price_large.jl

using Random
using Statistics
using Printf
using DiscreteMarket

const N_MARKETS = 500
const T_SIM     = 1_000
const MAX_UNITS = 5

# ── Surplus helper ────────────────────────────────────────────────────────────

function surplus_f64(m::GoodMarket, q::Int) :: Float64
    q <= 0 && return 0.0
    all_wtp  = sort([Float64(v) for d in m.consumers for v in d.wtp],  rev=true)
    all_wtac = sort([Float64(c) for f in m.firms     for c in f.wtac])
    n = min(q, length(all_wtp), length(all_wtac))
    n == 0 && return 0.0
    sum(all_wtp[1:n]) - sum(all_wtac[1:n])
end

# ── Precomputed BiGeo price structs (binary-search for speed) ─────────────────

struct BGIPConsumer
    log_p  :: Float64
    log_hi :: Float64
    wtp    :: Vector{Float64}   # active (≥ p*), sorted descending
end

struct BGIPFirm
    log_lo :: Float64
    log_p  :: Float64
    wtac   :: Vector{Float64}   # active (≤ p*), sorted ascending
end

function BGIPConsumer(d::ConsumerDemand, p_star::Price)
    lp     = log(Float64(p_star))
    active = sort([Float64(v) for v in d.wtp if v >= p_star], rev=true)
    lhi    = isempty(active) ? lp : log(active[1])
    BGIPConsumer(lp, lhi, active)
end

function BGIPFirm(f::FirmSupply, p_star::Price)
    lp     = log(Float64(p_star))
    active = sort([Float64(c) for c in f.wtac if c <= p_star])
    llo    = isempty(active) ? lp : log(active[1])
    BGIPFirm(llo, lp, active)
end

function bgi_demand(c::BGIPConsumer, rng::AbstractRNG) :: Int
    isempty(c.wtp) && return 0
    c.log_hi ≈ c.log_p && return length(c.wtp)
    pv = exp(c.log_p + rand(rng) * (c.log_hi - c.log_p))
    searchsortedlast(c.wtp, pv, rev=true)
end

function bgi_supply(f::BGIPFirm, rng::AbstractRNG) :: Int
    isempty(f.wtac) && return 0
    f.log_p ≈ f.log_lo && return length(f.wtac)
    pv = exp(f.log_lo + rand(rng) * (f.log_p - f.log_lo))
    searchsortedlast(f.wtac, pv)
end

# ── Record ────────────────────────────────────────────────────────────────────

struct CompLargeRecord
    q_wge        :: Int
    s_wge        :: Float64
    eff_std      :: Float64
    eff_bgip     :: Float64
    mean_q_std   :: Float64
    mean_q_bgip  :: Float64
end

# ── Per-size run ──────────────────────────────────────────────────────────────

function run_size(nc, nf; seed_base)
    records    = CompLargeRecord[]
    n_no_equil = 0

    for seed in 1:N_MARKETS
        rng = MersenneTwister(seed + seed_base)
        m, r = generate_good_market(rng;
            good=1, n_consumers=nc, n_firms=nf,
            max_units=MAX_UNITS, Q=100)

        r.cleared || (n_no_equil += 1; continue)
        p_star = r.price
        q_star = aggregate_demand(m, p_star)
        s_wge  = surplus_f64(m, q_star)
        s_wge > 0 || continue

        # ── Standard ZIT ─────────────────────────────────────────────────────
        zi_c = [ZIConsumer(1, p_star * length(d.wtp)) for d in m.consumers]
        zi_f = [ZIFirm(1, length(f.wtac))             for f in m.firms]
        zm   = ZIMarket(zi_c, zi_f, 1)
        res  = zi_simulate(zm, [p_star], T_SIM; seed=seed)
        t_std = res.traded[:, 1]
        eff_std = mean(surplus_f64(m, t) for t in t_std) / s_wge

        # ── BiGeo price-space ZIT ─────────────────────────────────────────────
        bgip_c = [BGIPConsumer(d, p_star) for d in m.consumers]
        bgip_f = [BGIPFirm(f, p_star)     for f in m.firms]
        rng2   = MersenneTwister(seed)
        t_bgip = Vector{Int}(undef, T_SIM)
        for t in 1:T_SIM
            d = sum(bgi_demand(c, rng2) for c in bgip_c; init=0)
            s = sum(bgi_supply(f, rng2) for f in bgip_f; init=0)
            t_bgip[t] = min(d, s)
        end
        eff_bgip = mean(surplus_f64(m, t) for t in t_bgip) / s_wge

        (isnan(eff_std) || isnan(eff_bgip)) && continue

        push!(records, CompLargeRecord(
            q_star, s_wge, eff_std, eff_bgip,
            mean(t_std), mean(Float64.(t_bgip))
        ))
    end

    records, n_no_equil
end

# ── Print ─────────────────────────────────────────────────────────────────────

function print_size(nc, nf, records, n_no_equil)
    n      = length(records)
    q_wge  = [r.q_wge       for r in records]
    s_wge  = [r.s_wge       for r in records]
    e_std  = [r.eff_std     for r in records]
    e_bgip = [r.eff_bgip    for r in records]
    q_std  = [r.mean_q_std  for r in records]
    q_bgip = [r.mean_q_bgip for r in records]

    println("=" ^ 62)
    @printf("%d consumers + %d firms  (%d agents total)\n", nc, nf, nc+nf)
    println("  max $MAX_UNITS units/agent, $N_MARKETS seeds, $T_SIM ZI trials each")
    println("=" ^ 62)
    @printf("Markets analyzed:  %d / %d  (%d had no WGE)\n\n", n, N_MARKETS, n_no_equil)

    println("── WGE outcomes ────────────────────────────────────────────")
    @printf("  q*:  mean=%.1f  median=%.1f  [%.0f, %.0f]\n",
            mean(q_wge), median(q_wge), minimum(q_wge), maximum(q_wge))
    @printf("  S*:  mean=%.1f  median=%.1f  [%.1f, %.1f]\n\n",
            mean(s_wge), median(s_wge), minimum(s_wge), maximum(s_wge))

    println("── Quantity means ──────────────────────────────────────────")
    @printf("  WGE q*:           mean=%.1f\n", mean(q_wge))
    @printf("  Standard ZIT q:   mean=%.1f  (%.1f%% of q*)\n",
            mean(q_std),  100*mean(q_std) /mean(q_wge))
    @printf("  BiGeo-Price q:    mean=%.1f  (%.1f%% of q*)\n\n",
            mean(q_bgip), 100*mean(q_bgip)/mean(q_wge))

    println("── Efficiency ──────────────────────────────────────────────")
    println("                Standard ZIT    BiGeo Price     Δ")
    @printf("  Mean:         %.4f          %.4f        %+.4f\n",
            mean(e_std), mean(e_bgip), mean(e_bgip)-mean(e_std))
    @printf("  Median:       %.4f          %.4f        %+.4f\n",
            median(e_std), median(e_bgip), median(e_bgip)-median(e_std))
    @printf("  Std:          %.4f          %.4f        %+.4f\n",
            std(e_std), std(e_bgip), std(e_bgip)-std(e_std))
    @printf("  Min:          %.4f          %.4f        %+.4f\n",
            minimum(e_std), minimum(e_bgip), minimum(e_bgip)-minimum(e_std))
    println()
    @printf("  Sharpe:       %.3f           %.3f\n\n",
            mean(e_std)/std(e_std), mean(e_bgip)/std(e_bgip))

    println("── Tail counts ─────────────────────────────────────────────")
    println("  Threshold │  Standard    BiGeo-Price")
    for thresh in [0.90, 0.95, 0.98, 0.99]
        ns = count(e -> e >= thresh, e_std)
        nb = count(e -> e >= thresh, e_bgip)
        @printf("  eff≥%.2f  │  %3d/%d (%.1f%%)  %3d/%d (%.1f%%)\n",
                thresh, ns, n, 100*ns/n, nb, n, 100*nb/n)
    end
    println()

    n_bgip_wins = count(i -> e_bgip[i] > e_std[i], 1:n)
    n_std_wins  = count(i -> e_bgip[i] < e_std[i], 1:n)
    println("── Head-to-head (same market, same seed) ───────────────────")
    @printf("  BiGeo-P > Std:  %d / %d  (%.1f%%)\n", n_bgip_wins, n, 100*n_bgip_wins/n)
    @printf("  Std > BiGeo-P:  %d / %d  (%.1f%%)\n\n", n_std_wins,  n, 100*n_std_wins/n)
end

# ── Run ───────────────────────────────────────────────────────────────────────

for (nc, nf, sb) in [(50, 50, 30_000), (250, 250, 31_000)]
    recs, no_eq = run_size(nc, nf; seed_base=sb)
    print_size(nc, nf, recs, no_eq)
end
