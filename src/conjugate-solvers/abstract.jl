abstract type AbstractConjugateSolver end

# Interface
set_odpairs(cmodel::AbstractConjugateSolver, commodities::AbstractVector{Commodity}) =
    set_odpairs(cmodel, [(comm.orig, comm.dest) for comm in commodities])

set_odpairs(cmodel::AbstractConjugateSolver, ks::AbstractVector{Int}) =
    set_odpairs(cmodel, cmodel.prob.K[ks])

function set_demands(cmodel::AbstractConjugateSolver, demands)
    set_demands(cmodel, demands, keys(demands))
end

# Calculate demands from paths
function calculate_demands(cmodel::AbstractConjugateSolver, paths)
    prob = problem(cmodel)
    nk, η = num_commodities(cmodel), weights(cmodel)
    a1 = tolled_arcs(prob)

    # Extract the set of tolled arcs for each path
    arcdict = srcdst_to_index(prob)
    a1set = BitSet(a1)
    tolled_sets = [BitSet(path_tolled_arcs(path, arcdict, a1set)) for path in paths]

    demands = Dict(a => sum(a in tolled_sets[k] ? η[k] : 0 for k in 1:nk) for a in a1)
    return demands
end

# Set paths
function set_paths(cmodel::AbstractConjugateSolver, paths, args...; set_odpairs=true)
    set_odpairs && NetPricing.set_odpairs(cmodel, [(path[1], path[end]) for path in paths])
    set_demands(cmodel, calculate_demands(cmodel, paths), args...)
end

# Bilevel feasible test
function is_bilevel_feasible(cmodel::AbstractConjugateSolver, paths; kwargs...)
    set_paths(cmodel, paths, BitSet(); kwargs...)   # The value of t is unimportant, so no priority for any arc
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

function solve(cmodel::AbstractConjugateSolver, demands, args...)
    set_demands(cmodel, demands, args...)
    return solve(cmodel)
end
