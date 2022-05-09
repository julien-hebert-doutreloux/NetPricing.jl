# Bilevel feasibility tester from conjugate solver
struct ConjugateBFTester{S<:AbstractConjugateSolver} <: AbstractBFTester
    solver::S
end

ConjugateBFTester{S}(prob::Problem, num_commodities::Int; kwargs...) where {S} = ConjugateBFTester(S(prob, num_commodities; kwargs...))
ConjugateBFTester(::Type{S}, prob::Problem, num_commodities::Int; kwargs...) where {S} = ConjugateBFTester{S}(prob, num_commodities; kwargs...)

problem(tester::ConjugateBFTester) = problem(tester.solver)
num_commodities(tester::ConjugateBFTester) = num_commodities(tester.solver)

# Set OD-pairs
set_odpairs(tester::ConjugateBFTester, odpairs::AbstractVector{Tuple{Int,Int}}) = set_odpairs(tester.solver, odpairs)

# Bilevel feasible test
function is_bilevel_feasible(tester::ConjugateBFTester, paths; set_odpairs=true)
    cmodel = tester.solver

    set_odpairs && NetPricing.set_odpairs(tester, [(path[1], path[end]) for path in paths])
    set_demands(cmodel, count_tolled_arcs(problem(tester), paths), BitSet())     # The value of t is unimportant, so no priority for any arc

    optimize!(cmodel)

    objval = objective_value(cmodel)
    arccosts = srcdst_to_cost(problem(cmodel))
    expecting = sum(get_path_cost(path, arccosts) * weight for (path, weight) in zip(paths, weights(cmodel)))

    return objval >= expecting - 1e-6
end
