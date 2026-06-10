# Phase 7 (optional extension): Padé approximation for irrational utility classes
#
# Load with:
#   include(joinpath(pkgdir(DiscreteMarket), "src", "extensions", "pade.jl"))
#
# Depends on: core.jl (Price, PriceVec, Bundle, RatVec)

# ── Abstract interface ─────────────────────────────────────────────────────────

abstract type RationalDemandSystem end

# Extend DiscreteMarket.walras_residual with new methods for RationalDemandSystem types.
# walras_residual is exported from DiscreteMarket; we import it so our new methods
# attach to the same generic function rather than shadowing it.
import DiscreteMarket: walras_residual

# demand(d, p, ω) :: Bundle — integer quantities given prices p and endowment ω
function demand end

# ── Cobb-Douglas (exact rational arithmetic) ───────────────────────────────────

struct CobbDouglas <: RationalDemandSystem
    α :: RatVec   # preference shares ∈ ℚ₊, must sum to 1
end

function demand(d::CobbDouglas, p::PriceVec, ω::Bundle) :: Bundle
    I = sum(p[j] * ω[j] for j in eachindex(p); init = 0//1)
    [floor(Int, d.α[j] * I / p[j]) for j in eachindex(p)]
end

function walras_residual(d::CobbDouglas, p::PriceVec, ω::Bundle) :: Rational{Int64}
    I = sum(p[j] * ω[j] for j in eachindex(p); init = 0//1)
    x = demand(d, p, ω)
    I - sum(p[j] * x[j] for j in eachindex(p); init = 0//1)
end

# ── Padé approximation utilities ──────────────────────────────────────────────

# Generalized binomial coefficient C(α, k) = α(α−1)⋯(α−k+1) / k!
function _gen_binom(α::Float64, k::Int) :: Float64
    k == 0 && return 1.0
    prod(α - i for i in 0:k-1) / factorial(k)
end

# Taylor coefficients of (1+t)^α at t = 0: result[k+1] = C(α, k)
function _taylor_power(α::Float64, n::Int) :: Vector{Float64}
    [_gen_binom(α, k) for k in 0:n]
end

# [m/m] Padé approximant from Taylor coefficients c = [c₀, c₁, ..., c_{2m}].
# Returns (p_coeffs, q_coeffs) where q_coeffs[1] = 1 (normalization).
# Condition: Q(t)·f(t) − P(t) = O(t^{2m+1}).
function pade_from_taylor(c::Vector{Float64}, m::Int) ::
        Tuple{Vector{Float64}, Vector{Float64}}
    length(c) >= 2m + 1 || error("need at least $(2m+1) Taylor coefficients, got $(length(c))")
    m == 0 && return ([c[1]], [1.0])

    # Hankel system: H[i,j] = c[m+i-j+1], rhs[i] = -c[m+i+1]
    H   = [c[m + i - j + 1] for i in 1:m, j in 1:m]
    rhs = [-c[m + i + 1]     for i in 1:m]
    q_tail = H \ rhs

    q = vcat(1.0, q_tail)
    p = [c[k+1] + sum(q[j+1] * c[k-j+1] for j in 1:k; init = 0.0) for k in 0:m]
    p, q
end

# ── 1D Padé approximant for p^α on a compact price domain ────────────────────
#
# Approximant: p^α  ≈  scale · P(t) / Q(t),   t = (p − center)/center
# where P, Q are polynomials of degree m with Q(0) = 1, rationalized to
# denominator bound Q_bound.

struct PadeApprox1D
    α      :: Float64
    m      :: Int
    center :: Float64    # midpoint of [p_min, p_max]
    scale  :: Float64    # center^α
    p_num  :: RatVec     # rationalized numerator coefficients
    p_den  :: RatVec     # rationalized denominator coefficients
end

function PadeApprox1D(α::Float64, m::Int,
                      p_min::Price, p_max::Price; Q_bound::Int = 1000)
    p0  = (Float64(p_min) + Float64(p_max)) / 2.0
    sc  = p0^α
    c   = _taylor_power(α, 2m)
    nf, df = pade_from_taylor(c, m)
    tol = 1.0 / Q_bound
    PadeApprox1D(α, m, p0, sc,
                 rationalize.(nf, tol = tol),
                 rationalize.(df, tol = tol))
end

# Float64 evaluation (used inside demand computations to avoid rational overflow)
function (pa::PadeApprox1D)(p::Price) :: Float64
    t   = (Float64(p) - pa.center) / pa.center
    num = evalpoly(t, Float64.(pa.p_num))
    den = evalpoly(t, Float64.(pa.p_den))
    pa.scale * num / den
end

# Exact rational evaluation — useful for error analysis; may overflow for large m
function eval_rational(pa::PadeApprox1D, p::Price) :: Rational{Int64}
    t   = (p - rationalize(pa.center, tol = 1.0 / 10^6)) //
          rationalize(pa.center,      tol = 1.0 / 10^6)
    num = sum(pa.p_num[k+1] * t^k for k in 0:pa.m)
    den = sum(pa.p_den[k+1] * t^k for k in 0:pa.m)
    rationalize(pa.scale, tol = 1.0 / 10^6) * num // den
end

# Exact Float64 value — reference for error measurement
exact_power(pa::PadeApprox1D, p::Price) :: Float64 = Float64(p)^pa.α

# ── CES demand with Padé approximation ────────────────────────────────────────
#
# Marshallian demand:  xⱼ = I · αⱼ^σ · pⱼ^(−σ) / Σₗ αₗ^σ · pₗ^(1−σ)
# Use the identity p^(1−σ) = p · p^(−σ), so only one PadeApprox1D is needed.
# This also avoids a degenerate Hankel matrix when 1−σ is a non-positive integer.

struct CESApprox <: RationalDemandSystem
    α        :: RatVec         # preference weights, Σα = 1
    α_pow    :: RatVec         # rationalized αⱼ^σ
    pade_neg :: PadeApprox1D   # approximates p^(−σ)
end

function CESApprox(α::RatVec, σ::Float64,
                   p_min::Price, p_max::Price;
                   m::Int = 2, Q_bound::Int = 1000)
    tol      = 1.0 / Q_bound
    α_pow    = rationalize.(Float64.(α) .^ σ, tol = tol)
    pade_neg = PadeApprox1D(-σ, m, p_min, p_max; Q_bound = Q_bound)
    CESApprox(α, α_pow, pade_neg)
end

function demand(d::CESApprox, p::PriceVec, ω::Bundle) :: Bundle
    I_f   = sum(Float64(p[j]) * ω[j] for j in eachindex(p))
    # denom = Σⱼ αⱼ^σ · pⱼ^(1−σ) = Σⱼ αⱼ^σ · pⱼ · pⱼ^(−σ)   (exact identity)
    denom = sum(Float64(d.α_pow[j]) * Float64(p[j]) * d.pade_neg(p[j])
                for j in eachindex(p))
    numer = [Float64(d.α_pow[j]) * d.pade_neg(p[j]) for j in eachindex(p)]
    [floor(Int, I_f * numer[j] / denom) for j in eachindex(p)]
end

function walras_residual(d::CESApprox, p::PriceVec, ω::Bundle) :: Rational{Int64}
    I = sum(p[j] * ω[j] for j in eachindex(p); init = 0//1)
    x = demand(d, p, ω)
    I - sum(p[j] * x[j] for j in eachindex(p); init = 0//1)
end

# Approximation error at a single price: |R(p) − p^α| / p^α
function relative_error(pa::PadeApprox1D, p::Price) :: Float64
    exact = exact_power(pa, p)
    abs(pa(p) - exact) / exact
end
