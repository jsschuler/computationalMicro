# Overnight test suite — large markets, high seed counts.
#
# Run with:
#   julia --project=. test/overnight.jl
#
# Estimated runtime: ~1.5–2 hours on a modern laptop (calibrated at 8 ms/XL market).
# T1–T3, T6: 100 k XL single-good markets each (~13 min each)
# T4:         30 k XL multi-good markets       (~16 min)
# T5:        100 k XL revealed-preference      (~20 min)
# T7:         50 k Phase 9 factor markets       (~17 min)
# The regular test suite (runtests.jl) runs in ~10 s and is unchanged.
#
# Parameter guide (adjust N_* constants to trade coverage for speed):
#   N_XL   = 100_000  →  ~15 min per single-good XL test set
#   N_MG   =  30_000  →  ~15 min per multi-good XL test set
#   N_P8   =  50_000  →   ~8 min for Phase 8 focal stats
#   N_P9   =  50_000  →  ~15 min for Phase 9 factor market stress

using Test
using Random
using Statistics
using DiscreteMarket

include(joinpath(pkgdir(DiscreteMarket), "src", "extensions", "factor_markets.jl"))
include(joinpath(pkgdir(DiscreteMarket), "src", "extensions", "stern_brocot.jl"))

# ── Scale constants ────────────────────────────────────────────────────────────

const N_XL = 100_000   # single-good XL tests
const N_MG =  30_000   # multi-good XL tests
const N_P8 =  50_000   # Phase 8 focal stats
const N_P9 =  50_000   # Phase 9 factor market

# ── Market size helpers ────────────────────────────────────────────────────────

# Extra-large single-good: 50–200 consumers, 50–150 firms, 10–30 units, Q=100_000
function xl_good_market(rng)
    generate_good_market(rng;
        good       = 1,
        n_consumers = rand(rng, 50:200),
        n_firms     = rand(rng, 50:150),
        max_units   = rand(rng, 10:30),
        Q           = 100_000)
end

# Extra-large multi-good: k=2:6 goods, 20–80 consumers, 20–60 firms, 5–20 units
function xl_market(seed)
    rng = MersenneTwister(seed)
    generate_market(seed;
        k           = rand(rng, 2:6),
        n_consumers = rand(rng, 20:80),
        n_firms     = rand(rng, 20:60),
        max_units   = rand(rng, 5:20))
end

# ── Test 1: find_equilibrium correctness on XL markets ───────────────────────

@testset "Overnight T1: find_equilibrium — XL markets ($N_XL seeds)" begin
    n_nonexist = 0
    for seed in 1:N_XL
        rng = MersenneTwister(seed)
        m, r = xl_good_market(rng)
        if r.cleared
            @test clears(m, r.price)
        else
            n_nonexist += 1
            @test !clears(m, r.price)
        end
    end
    frac = round(100 * n_nonexist / N_XL, digits=2)
    @info "T1 non-existence rate: $n_nonexist / $N_XL ($frac%)"
    @test n_nonexist < N_XL ÷ 5   # < 20% non-existence even in large markets
end

# ── Test 2: Walras' Law on XL single-good markets ────────────────────────────

@testset "Overnight T2: Walras' Law — XL single-good ($N_XL seeds)" begin
    n_tested = 0
    for seed in 1:N_XL
        rng = MersenneTwister(seed + 1_000_000)
        m, r = xl_good_market(rng)
        r.cleared || continue
        mkt = Market([m], 1)
        @test walras_residual(mkt, [r.price]) == 0//1
        n_tested += 1
    end
    @info "T2 Walras tested: $n_tested / $N_XL"
    @test n_tested >= N_XL ÷ 2
end

# ── Test 3: tatônnement ↔ find_equilibrium agreement on XL markets ───────────

@testset "Overnight T3: tatônnement agreement — XL markets ($N_XL seeds)" begin
    n_tested = 0
    for seed in 1:N_XL
        rng = MersenneTwister(seed + 2_000_000)
        m, r = xl_good_market(rng)
        r.cleared || continue
        p_tat = tatonnement(m)
        @test clears(m, r.price)
        @test clears(m, p_tat)
        n_tested += 1
    end
    @info "T3 tatônnement tested: $n_tested / $N_XL"
end

# ── Test 4: Walras' Law on XL multi-good markets ─────────────────────────────

@testset "Overnight T4: Walras' Law — XL multi-good ($N_MG seeds)" begin
    n_tested = 0
    for seed in 1:N_MG
        mkt, p_stars = xl_market(seed + 3_000_000)
        cleared = findall(r -> r.cleared, p_stars)
        isempty(cleared) && continue
        p_known  = Price[p_stars[j].price for j in cleared]
        residual = sum(p_known[i] * excess_demand(mkt.goods[cleared[i]], p_known[i])
                       for i in eachindex(cleared))
        @test residual == 0//1
        n_tested += 1
    end
    @info "T4 multi-good Walras tested: $n_tested / $N_MG"
    @test n_tested >= N_MG ÷ 2
end

# ── Test 5: revealed preference on XL markets ────────────────────────────────

@testset "Overnight T5: revealed preference — XL markets ($N_XL seeds)" begin
    for seed in 1:N_XL
        rng = MersenneTwister(seed + 4_000_000)
        m, r = xl_good_market(rng)
        r.cleared || continue
        for d in m.consumers
            @test check_revealed_preference(d, r.price)
        end
    end
end

# ── Test 6: Phase 8 focal stats on large markets ─────────────────────────────

@testset "Overnight T6: Phase 8 focal_stats — XL markets ($N_P8 seeds)" begin
    depth_savings = Int[]
    for seed in 1:N_P8
        rng = MersenneTwister(seed + 5_000_000)
        m, r = xl_good_market(rng)
        r.cleared || continue
        stats = focal_stats(m, 1//20)
        @test stats.depth_saving >= 0
        push!(depth_savings, stats.depth_saving)
    end
    @info "T6 focal stats: $(length(depth_savings)) markets; mean depth saving $(round(mean(depth_savings), digits=2))"
    @test length(depth_savings) >= N_P8 ÷ 2
end

# ── Test 7: Phase 9 factor market stress — large markets ─────────────────────

@testset "Overnight T7: Phase 9 factor market — large ($N_P9 seeds)" begin
    n_labor_cleared  = 0
    n_goods_nonexist = 0

    for seed in 1:N_P9
        rng = MersenneTwister(seed + 6_000_000)
        k = rand(rng, 1:5)

        goods = [begin
            nc = rand(rng, 2:20)
            consumers = [ConsumerDemand(j,
                sort([rand(rng, 1:50)//1 for _ in 1:rand(rng, 1:8)], rev=true))
                for _ in 1:nc]
            GoodMarket(j, consumers, FirmSupply[])
        end for j in 1:k]
        m = Market(goods, k)

        ℓs = [1//4, 1//2, 1//1, 3//2, 2//1, 3//1]
        firms_wl = [[FirmWithLabor(j, rand(rng, ℓs), rand(rng, 1:8))]
                    for j in 1:k]

        max_labor = sum(f.labor_per_unit * f.capacity
                        for j in 1:k for f in firms_wl[j]; init=0//1)
        endowment = max(1//4, max_labor * (rand(rng, 1:9) // 10))
        lm = LaborMarket(1, [endowment])

        p_stars, w_star = solve_wge_with_labor(m, firms_wl, lm)

        @test length(p_stars) == k
        @test w_star >= 1//100
        @test w_star <= 100//1

        for j in 1:k
            r = p_stars[j]
            r.cleared || (n_goods_nonexist += 1; continue)
            mc_j = marginal_cost(firms_wl[j][1], w_star)
            gm_j = GoodMarket(j, m.goods[j].consumers,
                              [FirmSupply(j, fill(mc_j, firms_wl[j][1].capacity))])
            @test clears(gm_j, r.price)
        end

        total_labor = sum(
            firms_wl[j][1].labor_per_unit *
                min(firms_wl[j][1].capacity,
                    aggregate_demand(
                        GoodMarket(j, m.goods[j].consumers,
                            [FirmSupply(j, fill(marginal_cost(firms_wl[j][1], w_star),
                                               firms_wl[j][1].capacity))]),
                        p_stars[j].price))
            for j in 1:k; init=0//1)
        total_labor == labor_supply(lm) && (n_labor_cleared += 1)
    end

    @info "T7 Phase 9: labor cleared $n_labor_cleared / $N_P9; goods non-exist $n_goods_nonexist"
    @test n_labor_cleared >= N_P9 ÷ 20
end
