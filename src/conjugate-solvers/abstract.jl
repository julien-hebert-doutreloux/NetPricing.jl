abstract type AbstractConjugateSolver end

# Interface
set_odpairs(cmodel::AbstractConjugateSolver, commodities::AbstractVector{Commodity}) =
    set_odpairs(cmodel, [(comm.orig, comm.dest) for comm in commodities])

set_odpairs(cmodel::AbstractConjugateSolver, ks::AbstractVector{Int}) =
    set_odpairs(cmodel, cmodel.prob.K[ks])

function set_demands(cmodel::AbstractConjugateSolver, demands, discount::Real = 1-1e-5)
    a1 = tolled_arcs(problem(cmodel))
    set_demands(cmodel, demands, Dict(a1 .=> discount))
end

# Set paths
function set_paths(cmodel::AbstractConjugateSolver, paths, discount = 1-1e-5; set_odpairs=true)
    prob = problem(cmodel)
    nk, η = num_commodities(cmodel), weights(cmodel)
    a1 = tolled_arcs(prob)

    # Set od pairs if required
    if set_odpairs
        NetPricing.set_odpairs(cmodel, [(path[1], path[end]) for path in paths])
    end

    # Extract the set of tolled arcs for each path
    arcdict = srcdst_to_index(prob)
    a1set = BitSet(a1)
    tolled_sets = [BitSet(path_tolled_arcs(path, arcdict, a1set)) for path in paths]

    # Set the demands
    demands = Dict(a => sum(a in tolled_sets[k] ? η[k] : 0 for k in 1:nk) for a in a1)
    set_demands(cmodel, demands, discount)
end

# Bilevel feasible test
function is_bilevel_feasible(cmodel::AbstractConjugateSolver, paths; kwargs...)
    set_paths(cmodel, paths, 1; kwargs...)  # The value of t is not important, hence no discount
    optimize!(cmodel)
    objval = objective_value(cmodel)
    arccosts = srcdst_to_cost(problem(cmodel))
    expecting = sum(get_path_cost(path, arccosts) * weight for (path, weight) in zip(paths, weights(cmodel)))
    return objval >= expecting - 1e-6
end

# Solve for g(w) and t
function solve(cmodel::AbstractConjugateSolver)
    optimize!(cmodel)
    return objective_value(cmodel), tvals(cmodel)
end

function solve(cmodel::AbstractConjugateSolver, demands, discount=1-1e-5)
    set_demands(cmodel, demands, discount)
    return solve(cmodel)
end
