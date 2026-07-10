# Decentralized ZIT — adaptive price beliefs, NO global signals.
#
# Agents know ONLY:
#   - Own WTP / WTAC array
#   - Own realized trade k_i from each period's pro-rata allocation
#   - Nothing else: p*, D, S, q*, or other agents' actions are NEVER observed
#
# Each agent maintains a private price belief p̂_i, initialized maximally:
#   Consumers: p̂_i = 0       → all own units look profitable (n_active = n_total)
#   Firms:     p̂_j = max_WTAC → all own units look profitable (n_active = n_total)
#
# Price belief updates (additive, step DELTA):
#   Consumer rationed (k < d):          p̂_i += DELTA   [demanded too many → price too low]
#   Consumer at ratchet cap, not ration: p̂_i -= DELTA  [may be too conservative → lower]
#   Consumer ratcheting up, not ration:  hold            [still climbing → gather info]
#   Firm rationed (k < d):              p̂_j -= DELTA   [supplied too many → price too high]
#   Firm at ratchet cap, not ration:    p̂_j += DELTA   [may be too conservative → raise]
#   Firm ratcheting up, not ration:     hold
#
# Quantity: ratchet capped at n_active(p̂_i) = count(WTP ≥ p̂_i)
#   d_i(t) = min(n_active_i(p̂_i(t)), k_i(t-1) + 1)
#
# Convergence mechanism:
#   Agents rationed repeatedly → raise/lower p̂ → n_active shrinks toward WGE optimum
#   At WGE: D = S, no rationing, at cap → minor p̂ oscillation of ±DELTA around p*
#   Efficiency loss from oscillation → 0 as DELTA → 0

using Random
using Statistics
using Printf
using DiscreteMarket

const DELTA   = 1.0    # additive price step per turn (units = price)
const T_TURNS = 300

# ── Market-level surplus ──────────────────────────────────────────────────────

function surplus_f64(m::GoodMarket, q::Int) :: Float64
    q <= 0 && return 0.0
    all_wtp  = sort([Float64(v) for d in m.consumers for v in d.wtp], rev=true)
    all_wtac = sort([Float64(c) for f in m.firms     for c in f.wtac])
    n = min(q, length(all_wtp), length(all_wtac))
    n == 0 && return 0.0
    sum(all_wtp[1:n]) - sum(all_wtac[1:n])
end

# ── Pro-rata realized trade ───────────────────────────────────────────────────

@inline realized_c(draw, D, S) = D == 0 ? 0 : (D <= S ? draw : floor(Int, draw * S / D))
@inline realized_f(draw, D, S) = S == 0 ? 0 : (S <= D ? draw : floor(Int, draw * D / S))

# ── Per-agent data (immutable; mutable state kept in separate arrays) ─────────

struct DecentConsumer
    wtp :: Vector{Float64}   # ALL valuations, sorted descending
end

struct DecentFirm
    wtac :: Vector{Float64}  # ALL costs, sorted ascending
end

DecentConsumer(d::ConsumerDemand) = DecentConsumer(sort(Float64.(d.wtp), rev=true))
DecentFirm(f::FirmSupply)         = DecentFirm(sort(Float64.(f.wtac)))

@inline n_active_c(c::DecentConsumer, p::Float64) = count(v -> v >= p, c.wtp)
@inline n_active_f(f::DecentFirm,     p::Float64) = count(v -> v <= p, f.wtac)

# ── Price belief update rules ─────────────────────────────────────────────────

@inline function update_p_c(p::Float64, d::Int, k::Int, na::Int) :: Float64
    k < d              && return p + DELTA          # rationed → price too low → raise
    na == 0            && return max(p - DELTA, 0.0) # all units priced out → lower to re-enter
    d == na && na > 0  && return max(p - DELTA, 0.0) # at ratchet cap, not rationed → maybe too high → lower
    return p                                          # ratcheting up, not rationed → hold
end

@inline function update_p_f(p::Float64, d::Int, k::Int, na::Int) :: Float64
    k < d              && return max(p - DELTA, 0.0) # rationed → price too high → lower
    na == 0            && return p + DELTA            # all units priced out → raise to re-enter
    d == na && na > 0  && return p + DELTA            # at ratchet cap, not rationed → maybe too low → raise
    return p                                          # ratcheting up, not rationed → hold
end

# ── Simulate ──────────────────────────────────────────────────────────────────

function simulate_decent(consumers::Vector{DecentConsumer},
                          firms::Vector{DecentFirm},
                          m::GoodMarket, s_wge::Float64, q_star::Int) :: Tuple{Vector{Float64},Vector{Float64}}

    nc, nf = length(consumers), length(firms)

    # Initialize price beliefs: maximum aggressiveness
    # Consumer p̂=0 → all units look profitable; Firm p̂=max_WTAC → all units look profitable
    p_hat_c  = zeros(Float64, nc)
    p_hat_f  = [firms[j].wtac[end] for j in 1:nf]   # wtac sorted asc → last = max
    prev_k_c = zeros(Int, nc)
    prev_k_f = zeros(Int, nf)

    effs   = Vector{Float64}(undef, T_TURNS)
    q_frac = Vector{Float64}(undef, T_TURNS)

    for t in 1:T_TURNS
        # Current active unit counts under private price beliefs
        na_c = [n_active_c(consumers[i], p_hat_c[i]) for i in 1:nc]
        na_f = [n_active_f(firms[j],     p_hat_f[j]) for j in 1:nf]

        # Ratchet demands: probe one more than last period, capped at n_active(p̂)
        ds_c = [min(na_c[i], prev_k_c[i] + 1) for i in 1:nc]
        ds_f = [min(na_f[j], prev_k_f[j] + 1) for j in 1:nf]

        D = sum(ds_c)
        S = sum(ds_f)
        q = min(D, S)

        effs[t]   = surplus_f64(m, q) / s_wge
        q_frac[t] = q_star > 0 ? q / q_star : 1.0

        ks_c = [realized_c(ds_c[i], D, S) for i in 1:nc]
        ks_f = [realized_f(ds_f[j], D, S) for j in 1:nf]

        for i in 1:nc
            p_hat_c[i]  = update_p_c(p_hat_c[i], ds_c[i], ks_c[i], na_c[i])
            prev_k_c[i] = ks_c[i]
        end
        for j in 1:nf
            p_hat_f[j]  = update_p_f(p_hat_f[j], ds_f[j], ks_f[j], na_f[j])
            prev_k_f[j] = ks_f[j]
        end
    end

    effs, q_frac
end

# ── Baseline: standard ZIT (no p* needed) ────────────────────────────────────

@inline function baseline_std(consumers, firms, m, s_wge, T, rng)
    nc, nf = length(consumers), length(firms)
    mean(begin
        dc = sum(rand(rng, 0:length(consumers[i].wtp)) for i in 1:nc; init=0)
        ds = sum(rand(rng, 0:length(firms[j].wtac))    for j in 1:nf; init=0)
        surplus_f64(m, min(dc, ds)) / s_wge
    end for _ in 1:T)
end

# ── Run across markets ────────────────────────────────────────────────────────

function run_decent(; nc, nf, n_markets, seed_base, label,
                     max_units=5, report_turns=[1,2,3,5,10,20,50,100,200,300])

    sum_eff   = zeros(T_TURNS)
    sum_eff2  = zeros(T_TURNS)
    sum_qfrac = zeros(T_TURNS)
    n_analyzed = 0
    std_effs   = Float64[]

    for seed in 1:n_markets
        rng = MersenneTwister(seed + seed_base)
        m, r = generate_good_market(rng;
            good=1, n_consumers=nc, n_firms=nf,
            max_units=max_units, Q=100)
        r.cleared || continue
        p_star = r.price
        q_star = aggregate_demand(m, p_star)
        s_wge  = surplus_f64(m, q_star)
        s_wge > 0 || continue

        consumers = [DecentConsumer(d) for d in m.consumers]
        firms     = [DecentFirm(f)     for f in m.firms]

        effs, qf = simulate_decent(consumers, firms, m, s_wge, q_star)
        for t in 1:T_TURNS
            sum_eff[t]   += effs[t]
            sum_eff2[t]  += effs[t]^2
            sum_qfrac[t] += qf[t]
        end
        n_analyzed += 1

        rng2 = MersenneTwister(seed + seed_base + 99_000)
        push!(std_effs, baseline_std(consumers, firms, m, s_wge, T_TURNS, rng2))
    end

    n = n_analyzed
    println("=" ^ 70)
    println(label)
    @printf("  n=%d  T=%d  DELTA=%.1f\n\n", n, T_TURNS, DELTA)
    @printf("  Baseline — Static Std:   mean eff = %.4f\n",   mean(std_effs))
    @printf("  Target   — WGE:          mean eff = 1.0000\n")
    println()
    println("  Turn │  q/q*    Eff_mean  Eff_std  vs Std   vs WGE")
    println("  ─────┼──────────────────────────────────────────────")
    for t in report_turns
        t > T_TURNS && continue
        eff = sum_eff[t] / n
        qf  = sum_qfrac[t] / n
        sd  = sqrt(max(0.0, sum_eff2[t]/n - eff^2))
        @printf("  t=%3d │  %.4f  %.4f    %.4f   %+.4f  %+.4f\n",
                t, qf, eff, sd, eff - mean(std_effs), eff - 1.0)
    end
    println()
    final_eff = sum_eff[T_TURNS]  / n
    final_sd  = sqrt(max(0.0, sum_eff2[T_TURNS]/n - final_eff^2))
    final_qf  = sum_qfrac[T_TURNS] / n
    @printf("  Final (t=%d):\n", T_TURNS)
    @printf("    q/q* = %.4f\n", final_qf)
    @printf("    Efficiency: mean=%.4f  std=%.4f  vs WGE=%+.4f\n",
            final_eff, final_sd, final_eff - 1.0)
    println()

    (sum_eff ./ n, sum_qfrac ./ n, n)
end

# ── Main ──────────────────────────────────────────────────────────────────────

println("DECENTRALIZED ZIT — adaptive price beliefs, no global signals")
println("  Agents observe ONLY: own WTP/WTAC + own realized trade k_i per period")
println("  p*, D, S, q* are NEVER observed")
println()
println("  Price belief update (additive, step DELTA=$(DELTA)):")
println("    Consumer rationed:         p̂ += DELTA  (demanded too many → price too low)")
println("    Consumer at cap, clear:    p̂ -= DELTA  (may be too conservative → lower)")
println("    Consumer ratcheting:       hold         (still exploring, no signal)")
println("    Firm rationed:             p̂ -= DELTA  (supplied too many → price too high)")
println("    Firm at cap, clear:        p̂ += DELTA  (may be too conservative → raise)")
println("  Quantity:  d_i(t) = min(count(WTP ≥ p̂_i), k_i(t-1) + 1)  [ratchet]")
println("  Init:      consumer p̂=0 (all units active), firm p̂=max_WTAC (all units active)")
println()

RT = [1, 2, 3, 5, 10, 20, 50, 100, 200, 300]

# Small / mixed markets
let
    SEED_BASE  = 40_000
    sum_eff    = zeros(T_TURNS); sum_eff2  = zeros(T_TURNS)
    sum_qfrac  = zeros(T_TURNS)
    n          = 0
    std_effs   = Float64[]

    for seed in 1:1_000
        rng  = MersenneTwister(seed + SEED_BASE)
        n_c  = rand(rng, 1:8); n_f = rand(rng, 1:8); max_u = rand(rng, 1:6)
        m, r = generate_good_market(rng;
            good=1, n_consumers=n_c, n_firms=n_f, max_units=max_u, Q=100)
        r.cleared || continue
        p_star = r.price
        q_star = aggregate_demand(m, p_star)
        s_wge  = surplus_f64(m, q_star)
        s_wge > 0 || continue

        consumers = [DecentConsumer(d) for d in m.consumers]
        firms     = [DecentFirm(f)     for f in m.firms]

        effs, qf = simulate_decent(consumers, firms, m, s_wge, q_star)
        for t in 1:T_TURNS
            sum_eff[t]  += effs[t]; sum_eff2[t]  += effs[t]^2
            sum_qfrac[t] += qf[t]
        end
        n += 1

        rng2 = MersenneTwister(seed + SEED_BASE + 99_000)
        push!(std_effs, baseline_std(consumers, firms, m, s_wge, T_TURNS, rng2))
    end

    println("=" ^ 70)
    println("Small/mixed markets  (1–8C, 1–8F, 1–6 max_units)  n=$n / 1000")
    @printf("  T=%d  DELTA=%.1f\n\n", T_TURNS, DELTA)
    @printf("  Baseline — Static Std:   mean eff = %.4f\n", mean(std_effs))
    println()
    println("  Turn │  q/q*    Eff_mean  Eff_std  vs Std   vs WGE")
    println("  ─────┼──────────────────────────────────────────────")
    for t in RT
        t > T_TURNS && continue
        eff = sum_eff[t] / n; qf = sum_qfrac[t] / n
        sd  = sqrt(max(0.0, sum_eff2[t]/n - eff^2))
        @printf("  t=%3d │  %.4f  %.4f    %.4f   %+.4f  %+.4f\n",
                t, qf, eff, sd, eff - mean(std_effs), eff - 1.0)
    end
    println()
    fe = sum_eff[T_TURNS]/n; fq = sum_qfrac[T_TURNS]/n
    @printf("  Final (t=%d): q/q*=%.4f  eff=%.4f  std=%.4f  vs WGE=%+.4f\n\n",
            T_TURNS, fq, fe, sqrt(max(0.0,sum_eff2[T_TURNS]/n - fe^2)), fe - 1.0)
end

run_decent(nc=50, nf=50, n_markets=300, seed_base=50_000,
           label="100 agents (50C + 50F)  n=300 markets",
           report_turns=RT)

run_decent(nc=250, nf=250, n_markets=150, seed_base=51_000,
           label="500 agents (250C + 250F)  n=150 markets",
           report_turns=RT)
