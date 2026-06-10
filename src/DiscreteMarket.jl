module DiscreteMarket

using Random
using Statistics

# ── Phase 1: core types and aggregate functions ───────────────────────────────
include("core.jl")
export Price, PriceVec, Bundle, RatVec
export ConsumerDemand, FirmSupply, GoodMarket, Market
export supply_correspondence, aggregate_demand, aggregate_supply
export clears, excess_demand, walras_residual
export utility, marginal_utility, check_revealed_preference

# ── Phase 3: equilibrium solvers ─────────────────────────────────────────────
include("solvers.jl")
export find_equilibrium, tatonnement, solve_wge_exact, solve_wge

# ── Phase 2: random market generation ────────────────────────────────────────
include("generation.jl")
export draw_wtp, draw_wtac, generate_good_market, generate_market

# ── Phase 4: Zero Intelligence traders ───────────────────────────────────────
include("zi.jl")
export ZIConsumer, ZIFirm, ZIMarket
export zi_demand, zi_demand_mean, zi_supply, zi_supply_mean
export zi_period, zi_simulate

# ── Phases 5–6: comparison statistics and multi-market summary ───────────────
include("comparison.jl")
export total_surplus, wge_surplus, zi_efficiency, zi_clearing_price
export ComparisonStats, compare
export MultiMarketStats, multi_compare

end # module DiscreteMarket
