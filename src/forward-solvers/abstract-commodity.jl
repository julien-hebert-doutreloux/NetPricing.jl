abstract type CommodityForwardSolver end

# Empty solver
mutable struct CommodityForwardEmptySolver <: CommodityForwardSolver
    prob::AbstractCommodityProblem
end

problem(fmodel::CommodityForwardEmptySolver) = fmodel.prob
set_tolls(::CommodityForwardEmptySolver, tolls) = nothing
JuMP.optimize!(::CommodityForwardEmptySolver) = nothing
JuMP.objective_value(::CommodityForwardEmptySolver) = 0.
optimal_tolled_set(::CommodityForwardEmptySolver) = BitSet()
