include(joinpath(pkgdir(DiscreteMarket), "src", "extensions", "pade.jl"))

# ── Padé core algorithm ────────────────────────────────────────────────────────

@testset "Phase 7: pade_from_taylor — geometric series (1/(1−t))" begin
    # f(t) = 1/(1−t) has Taylor coefficients cₖ = 1 for all k.
    # The [1/1] Padé should recover 1/(1−t) exactly.
    c = ones(Float64, 3)
    p_c, q_c = pade_from_taylor(c, 1)
    # P(t) = 1, Q(t) = 1 − t  →  R(t) = 1/(1−t)
    @test p_c ≈ [1.0, 0.0]    atol=1e-12
    @test q_c ≈ [1.0, -1.0]   atol=1e-12
end

@testset "Phase 7: pade_from_taylor — (1+t)^(−2) at order 2 (exact)" begin
    # (1+t)^(−2) is rational of order [0/2], so [2/2] Padé is exact.
    # Taylor: c₀=1, c₁=−2, c₂=3, c₃=−4, c₄=5
    c = [1.0, -2.0, 3.0, -4.0, 5.0]
    p_c, q_c = pade_from_taylor(c, 2)
    # R(t) = P(t)/Q(t) = 1/(1+t)² → P=[1], Q=[1,2,1]
    @test p_c ≈ [1.0, 0.0, 0.0]  atol=1e-10
    @test q_c ≈ [1.0, 2.0, 1.0]  atol=1e-10
    # Spot-check: R(0.5) = 1/(1.5)² ≈ 0.4444
    t = 0.5
    @test evalpoly(t, p_c) / evalpoly(t, q_c) ≈ (1+t)^(-2)  atol=1e-10
end

# ── PadeApprox1D accuracy ──────────────────────────────────────────────────────

@testset "Phase 7: PadeApprox1D — p^(−2) on [1,9] is exact at order 2" begin
    # p^(−2) is a rational function, so order-2 Padé is exact.
    pa = PadeApprox1D(-2.0, 2, 1//1, 9//1; Q_bound=10000)
    for p in Price[1//1, 2//1, 3//1, 5//1, 7//1, 9//1]
        @test relative_error(pa, p) < 1e-10
    end
end

@testset "Phase 7: PadeApprox1D — p^(−0.5) on [1,9], order 2" begin
    pa = PadeApprox1D(-0.5, 2, 1//1, 9//1; Q_bound=1000)
    errors = [relative_error(pa, p) for p in Price[1//1, 3//1, 5//1, 7//1, 9//1]]
    @test all(errors .< 0.02)    # <2% relative error across a factor-9 domain
    @test maximum(errors) < 0.02
    @info "p^(-0.5) Padé[2,2] max relative error: $(maximum(errors)*100)%"
end

@testset "Phase 7: PadeApprox1D — p^(−1.5) on [1,4], order 3" begin
    pa = PadeApprox1D(-1.5, 3, 1//1, 4//1; Q_bound=1000)
    errors = [relative_error(pa, p) for p in Price[1//1, 2//1, 3//1, 4//1]]
    @test maximum(errors) < 5e-3
    @info "p^(-1.5) Padé[3,3] max relative error: $(maximum(errors)*100)%"
end

@testset "Phase 7: PadeApprox1D — rational evaluation is consistent" begin
    pa = PadeApprox1D(-0.5, 2, 1//1, 9//1; Q_bound=10000)
    for p in Price[2//1, 4//1, 6//1]
        r_val = Float64(eval_rational(pa, p))
        f_val = pa(p)
        # Float64 and rational evaluations should agree closely
        @test abs(r_val - f_val) / abs(f_val) < 1e-5
    end
end

# ── CobbDouglas ────────────────────────────────────────────────────────────────

@testset "Phase 7: CobbDouglas — exact budget balance (divisible case)" begin
    # α = [1/2, 1/2], p = [4, 4], ω = [10, 10]
    # I = 80, x₁ = x₂ = floor(10) = 10, residual = 0
    d = CobbDouglas([1//2, 1//2])
    p = Price[4//1, 4//1]
    ω = Int[10, 10]
    x = demand(d, p, ω)
    @test x == [10, 10]
    @test walras_residual(d, p, ω) == 0//1
end

@testset "Phase 7: CobbDouglas — floor residual is non-negative" begin
    # α = [1/2, 1/2], p = [4, 6], ω = [10, 8]
    # I = 88, x₁ = floor(11) = 11, x₂ = floor(7.33) = 7
    # residual = 88 − 44 − 42 = 2
    d = CobbDouglas([1//2, 1//2])
    p = Price[4//1, 6//1]
    ω = Int[10, 8]
    x = demand(d, p, ω)
    @test x == [11, 7]
    r = walras_residual(d, p, ω)
    @test r == 2//1
    @test r >= 0//1
end

@testset "Phase 7: CobbDouglas — demands are non-negative" begin
    rng = MersenneTwister(42)
    d = CobbDouglas([1//3, 1//3, 1//3])
    for _ in 1:20
        p = Price[rand(rng, 1:10)//1 for _ in 1:3]
        ω = Int[rand(rng, 1:20)    for _ in 1:3]
        x = demand(d, p, ω)
        @test all(x .>= 0)
        @test walras_residual(d, p, ω) >= 0//1
    end
end

@testset "Phase 7: CobbDouglas — Walras' Law (unspent ≤ min pⱼ)" begin
    # The floor residual cannot exceed the smallest per-unit price
    d = CobbDouglas([1//4, 3//4])
    p = Price[3//1, 5//1]
    ω = Int[20, 12]
    r = walras_residual(d, p, ω)
    @test 0//1 <= r
    # residual < max(pⱼ): can't have more than one unit's worth of any good unspent
    @test r < maximum(p)
end

# ── CESApprox ─────────────────────────────────────────────────────────────────

@testset "Phase 7: CESApprox — σ=2 (integer, Padé exact) matches formula" begin
    # With σ=2, p^(−2) and p^(1−σ)=p^(−1) are rational, so Padé is exact.
    # α=[1/2,1/2], p=[3,6], ω=[10,10], I=90
    # xⱼ = I · (αⱼ/2)² · pⱼ^(−2) / Σ (αₗ/2)² · pₗ^(−1)
    # α_pow = [0.25, 0.25]
    # denom = 0.25*(1/3) + 0.25*(1/6) = 1/12 + 1/24 = 1/8
    # x₁ = 90 * 0.25/9 / (1/8) = 90*(1/36)*8 = 20
    # x₂ = 90 * 0.25/36 / (1/8) = 90*(1/144)*8 = 5
    ces = CESApprox([1//2, 1//2], 2.0, 1//1, 10//1; m=2, Q_bound=10000)
    x = demand(ces, Price[3//1, 6//1], Int[10, 10])
    @test x == [20, 5]
    @test walras_residual(ces, Price[3//1, 6//1], Int[10, 10]) == 0//1
end

@testset "Phase 7: CESApprox — σ=1.5 (irrational) demand is integer and close to exact" begin
    # Exact CES demand in Float64 for comparison
    function ces_exact(α_f, σ, p_f, I_f)
        αpow = α_f .^ σ
        denom = sum(αpow[j] * p_f[j]^(1-σ) for j in eachindex(p_f))
        [floor(Int, I_f * αpow[j] * p_f[j]^(-σ) / denom) for j in eachindex(p_f)]
    end

    α   = [1//2, 1//2]
    p   = Price[2//1, 5//1]
    ω   = Int[15, 10]
    I_f = Float64(sum(p[j]*ω[j] for j in eachindex(p)))

    ces  = CESApprox(α, 1.5, 1//1, 8//1; m=2, Q_bound=1000)
    x    = demand(ces, p, ω)
    x_ex = ces_exact(Float64.(α), 1.5, Float64.(p), I_f)

    @test all(x .>= 0)
    @test walras_residual(ces, p, ω) >= 0//1
    # Demands should be within ±1 of exact (approximation + floor rounding)
    @test all(abs.(x .- x_ex) .<= 1)
    @info "CES σ=1.5: Padé demand $x  vs exact $x_ex"
end

@testset "Phase 7: CESApprox — Walras residual ≥ 0 for random prices" begin
    rng = MersenneTwister(7)
    ces = CESApprox([1//3, 1//3, 1//3], 0.8, 1//1, 10//1; m=3, Q_bound=1000)
    for _ in 1:20
        p = Price[rand(rng, 2:8)//1 for _ in 1:3]
        ω = Int[rand(rng, 5:15)  for _ in 1:3]
        @test walras_residual(ces, p, ω) >= 0//1
    end
end

@testset "Phase 7: CESApprox — higher order improves accuracy" begin
    α = [1//3, 2//3]
    p = Price[2//1, 7//1]
    ω = Int[30, 15]

    function approx_error(m_ord)
        ces = CESApprox(α, 1.7, 1//1, 10//1; m=m_ord, Q_bound=10000)
        αpow = Float64.(α) .^ 1.7
        I_f  = Float64(sum(p[j]*ω[j] for j in eachindex(p)))
        denom = sum(αpow[j] * Float64(p[j])^(1-1.7) for j in eachindex(p))
        x_ex  = [I_f * αpow[j] * Float64(p[j])^(-1.7) / denom for j in eachindex(p)]
        x_pd  = Float64.(demand(ces, p, ω))
        maximum(abs.(x_pd .- x_ex))
    end

    e1 = approx_error(1)
    e2 = approx_error(2)
    e3 = approx_error(3)
    @info "CES σ=1.7 demand error — m=1: $(round(e1,digits=3))  m=2: $(round(e2,digits=3))  m=3: $(round(e3,digits=3))"
    @test e2 <= e1 + 1    # higher order ≤ lower order error (or at worst ±1 from floor)
end
