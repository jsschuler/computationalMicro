# Phase 5: comparison statistics
# Depends on: core.jl

function total_surplus(m::GoodMarket, q::Int) :: Rational{Int64}
    q <= 0 && return 0//1
    # Efficient allocation: assign units to highest-WTP consumers and lowest-WTA firms
    all_wtp  = sort(vcat([d.wtp  for d in m.consumers]...), rev=true)
    all_wtac = sort(vcat([f.wtac for f in m.firms    ]...))
    n = min(q, length(all_wtp), length(all_wtac))
    n == 0 && return 0//1
    sum(all_wtp[1:n]) - sum(all_wtac[1:n])
end

function wge_surplus(m::GoodMarket, p_star::Price) :: Rational{Int64}
    q_star = aggregate_demand(m, p_star)
    total_surplus(m, q_star)
end

function zi_efficiency(m::GoodMarket, p_star::Price,
                       traded::Vector{Int}) :: Float64
    s_wge = Float64(wge_surplus(m, p_star))
    s_zi  = mean(Float64(total_surplus(m, q)) for q in traded)
    s_wge > 0 ? s_zi / s_wge : NaN
end

function zi_clearing_price(m::GoodMarket, zi_supply::Int) :: Union{Price, Nothing}
    candidates = sort(unique(vcat([d.wtp for d in m.consumers]...)))
    for p in candidates
        aggregate_demand(m, p) == zi_supply && return p
    end
    nothing
end

struct ComparisonStats
    good            :: Int
    p_wge           :: Price
    q_wge           :: Int
    zi_mean_demand  :: Float64
    zi_mean_supply  :: Float64
    zi_efficiency   :: Float64
    n_trials        :: Int
end

function compare(m::GoodMarket, p_wge::Price,
                 zi_d::Vector{Int}, zi_s::Vector{Int}) :: ComparisonStats
    ComparisonStats(
        m.good, p_wge,
        aggregate_demand(m, p_wge),
        mean(zi_d), mean(zi_s),
        zi_efficiency(m, p_wge, min.(zi_d, zi_s)),
        length(zi_d)
    )
end

# ── Phase 6: multi-good comparison ───────────────────────────────────────────

struct MultiMarketStats
    k           :: Int
    p_wge       :: Vector{Union{Price,Nothing}}
    q_wge       :: Vector{Int}
    surplus_wge :: Vector{Rational{Int64}}
    walras      :: Rational{Int64}    # p · Z(p*) — should be 0
end

function multi_compare(m::Market,
                       p_wge::Vector{Union{Price,Nothing}}) :: MultiMarketStats
    q_wge = Int[]
    surplus = Rational{Int64}[]
    p_known = Price[]   # prices where equilibrium exists, for Walras check

    for j in 1:m.k
        p = p_wge[j]
        if isnothing(p)
            push!(q_wge, 0)
            push!(surplus, 0//1)
        else
            push!(q_wge, aggregate_demand(m.goods[j], p))
            push!(surplus, wge_surplus(m.goods[j], p))
            push!(p_known, p)
        end
    end

    # Walras' Law on the subset of goods with known equilibrium prices
    known_idx = findall(!isnothing, p_wge)
    w = isempty(known_idx) ? 0//1 :
        sum(p_wge[j] * excess_demand(m.goods[j], p_wge[j])
            for j in known_idx)

    MultiMarketStats(m.k, p_wge, q_wge, surplus, w)
end
