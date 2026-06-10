using Test
using Random
using DiscreteMarket

# Run all phases, or restrict with: PHASES=1 julia --project=. test/runtests.jl
phases = let s = get(ENV, "PHASES", "")
    isempty(s) ? Set(1:8) : Set(parse.(Int, split(s, ",")))
end

1 in phases && include("phase1.jl")
2 in phases && include("phase2.jl")
3 in phases && include("phase3.jl")
4 in phases && include("phase4.jl")
(5 in phases || 6 in phases) && include("phase5.jl")
6 in phases && include("randomized.jl")
7 in phases && include("phase7.jl")
8 in phases && include("phase8.jl")
