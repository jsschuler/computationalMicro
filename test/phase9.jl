include(joinpath(pkgdir(DiscreteMarket), "src", "extensions", "factor_markets.jl"))

# ── Unit tests for helper functions ──────────────────────────────────────────

@testset "Phase 9: marginal_cost and supply_correspondence" begin
    f = FirmWithLabor(1, 1//2, 4)   # good 1, ℓ=1/2, capacity=4

    # At wage w=4, mc = 4 * 1/2 = 2
    @test marginal_cost(f, 4//1) == 2//1

    # p=3 > mc=2: firm strictly profitable → supply at capacity
    @test supply_correspondence(f, 3//1, 4//1) == (4, 4)
    # p=2 = mc=2: firm indifferent → (0, capacity)
    @test supply_correspondence(f, 2//1, 4//1) == (0, 4)
    # p=1 < mc=2: firm unprofitable → (0, 0)
    @test supply_correspondence(f, 1//1, 4//1) == (0, 0)
end

@testset "Phase 9: labor_supply and income" begin
    lm = LaborMarket(2, [1//2, 3//2])
    @test labor_supply(lm) == 2//1

    @test income(3//4, 8//1) == 6//1
end

@testset "Phase 9: labor_demand (single-good helper)" begin
    f1 = FirmWithLabor(1, 1//1, 3)
    f2 = FirmWithLabor(1, 2//1, 2)

    # At p=5, w=4: mc1=4 < p → supplies at capacity 3; mc2=8 > p → supplies 0
    @test labor_demand([f1, f2], 5//1, 4//1) == 3//1

    # At p=10, w=4: mc1=4 < p, mc2=8 < p → both at capacity; labor=1*3+2*2=7
    @test labor_demand([f1, f2], 10//1, 4//1) == 7//1

    # At p=4, w=4: mc1=4 = p → (0,3); mc2=8 > p → (0,0)
    # hi of f1 = 3, hi of f2 = 0 → total 3
    @test labor_demand([f1, f2], 4//1, 4//1) == 3//1
end

# ── End-to-end: single-good factor market ────────────────────────────────────
#
# Setup: consumer with WTP=[10,8,6,4], one labor firm ℓ=1 cap=4, endowment=2.
# Analysis:
#   At wage w, mc=w. Goods market: find_equilibrium picks p*=w where D(w)=q*.
#   D(w) = 2 when w ∈ (6, 8].  Labor = 1 * q* = q*.  Labor = 2 = endowment.
# The Stern-Brocot path from [1/100, 100] converges to some w ∈ (6,8] exactly.

@testset "Phase 9: single-good labor equilibrium" begin
    m1 = Market([
        GoodMarket(1, [ConsumerDemand(1, [10//1, 8//1, 6//1, 4//1])], FirmSupply[])
    ], 1)
    firms_wl = [[FirmWithLabor(1, 1//1, 4)]]
    lm1 = LaborMarket(1, [2//1])

    p_stars, w_star = solve_wge_with_labor(m1, firms_wl, lm1)

    @test length(p_stars) == 1
    @test p_stars[1].cleared

    # Reconstruct the goods market at the returned wage and verify clearing
    mc = marginal_cost(FirmWithLabor(1, 1//1, 4), w_star)
    gm = GoodMarket(1, m1.goods[1].consumers, [FirmSupply(1, fill(mc, 4))])
    @test clears(gm, p_stars[1].price)

    # Labor market must clear exactly
    q_star = aggregate_demand(gm, p_stars[1].price)
    @test 1//1 * q_star == labor_supply(lm1)
end

# ── End-to-end: two-good factor market ───────────────────────────────────────
#
# Setup:
#   Good 1: WTP=[10,6,2], labor firm ℓ=1, cap=3
#   Good 2: WTP=[8,4],    labor firm ℓ=2, cap=2
#   Endowment=4 (one worker supplies 4 units of labor).
#
# Analysis: at w ∈ (2, 4]:
#   mc1 = w ∈ (2,4] → D1(mc1) = 2 (10,6 ≥ w)
#   mc2 = 2w ∈ (4,8] → D2(mc2) = 1 (8 ≥ 2w for w ≤ 4)
#   Labor = 1*2 + 2*1 = 4 = endowment  ✓

@testset "Phase 9: two-good labor equilibrium" begin
    m2 = Market([
        GoodMarket(1, [ConsumerDemand(1, [10//1, 6//1, 2//1])], FirmSupply[]),
        GoodMarket(2, [ConsumerDemand(2, [8//1, 4//1])],         FirmSupply[])
    ], 2)
    firms_wl2 = [
        [FirmWithLabor(1, 1//1, 3)],
        [FirmWithLabor(2, 2//1, 2)]
    ]
    lm2 = LaborMarket(1, [4//1])

    p_stars2, w_star2 = solve_wge_with_labor(m2, firms_wl2, lm2)

    @test length(p_stars2) == 2
    @test p_stars2[1].cleared
    @test p_stars2[2].cleared

    # Both goods markets clear at returned prices
    mc1 = marginal_cost(FirmWithLabor(1, 1//1, 3), w_star2)
    mc2 = marginal_cost(FirmWithLabor(2, 2//1, 2), w_star2)
    gm1 = GoodMarket(1, m2.goods[1].consumers, [FirmSupply(1, fill(mc1, 3))])
    gm2 = GoodMarket(2, m2.goods[2].consumers, [FirmSupply(2, fill(mc2, 2))])
    @test clears(gm1, p_stars2[1].price)
    @test clears(gm2, p_stars2[2].price)

    # Labor market clears
    q1 = aggregate_demand(gm1, p_stars2[1].price)
    q2 = aggregate_demand(gm2, p_stars2[2].price)
    @test 1//1 * q1 + 2//1 * q2 == labor_supply(lm2)
end

# ── Mixed market: regular firms + labor firms ────────────────────────────────
#
# Some firms have exogenous costs (FirmSupply), some have labor-dependent costs.
# Regular firm: WTA=[2,2] (capacity 2 at cost 2).
# Labor firm: ℓ=1, cap=3.
# Consumer WTP=[10,8,6,4,2,1].  Endowment=3.
#
# Regular supply is always 2 units at p ≥ 2 (exogenous).
# Labor supply = capacity if p > w, 0 if p < w.
# At w=4: mc=4. D(4)=4 (10,8,6,4≥4). Regular supply=2, labor supply=2 (p=mc).
#   D=4 ∈ [S_reg_lo + S_lab_lo, S_reg_hi + S_lab_hi] = [2+0, 2+3] = [2,5].  Clears.
#   q=4. Labor from labor firm = 1*min(3,4-2)=1*2=2. Wait, q from labor firm=4-2=2.
#   Hmm, let me recalculate properly.
#   aggregate_demand = 4. Regular S = (2,2). Labor S = (0,3) at p=mc.
#   Total S = (2, 5). D=4 ∈ [2,5]. Clears.
#   Labor used by labor firm: min(cap=3, q=4) where q is total traded.
#   But regular firm supplies 2, so labor firm supplies 2 (to reach total 4).
#
# To keep the test simple, I verify the economic invariants rather than exact quantities.

@testset "Phase 9: mixed market (regular + labor firms)" begin
    m3 = Market([
        GoodMarket(1,
            [ConsumerDemand(1, [10//1, 8//1, 6//1, 4//1, 2//1, 1//1])],
            [FirmSupply(1, [2//1, 2//1])])    # 2 units at cost 2
    ], 1)
    firms_wl3 = [[FirmWithLabor(1, 1//1, 3)]]
    lm3 = LaborMarket(1, [3//1])

    p_stars3, w_star3 = solve_wge_with_labor(m3, firms_wl3, lm3)

    @test length(p_stars3) == 1
    @test p_stars3[1].cleared

    # Reconstruct and verify goods clearing
    mc3 = marginal_cost(FirmWithLabor(1, 1//1, 3), w_star3)
    gm3 = GoodMarket(1, m3.goods[1].consumers,
                     vcat(m3.goods[1].firms, [FirmSupply(1, fill(mc3, 3))]))
    @test clears(gm3, p_stars3[1].price)

    # Labor market clears: residual after regular firm supply equals endowment
    q3    = aggregate_demand(gm3, p_stars3[1].price)
    s_reg3 = sum(supply_correspondence(f, p_stars3[1].price)[2]
                 for f in m3.goods[1].firms; init=0)
    q_labor3 = max(0, q3 - s_reg3)
    @test 1//1 * min(3, q_labor3) == labor_supply(lm3)
end
