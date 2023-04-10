abstract type AbstractConjugateSolver end

# Interface
set_odpairs(cmodel::AbstractConjugateSolver, commodities::AbstractVector{Commodity}) =
    set_odpairs(cmodel, [(comm.orig, comm.dest) for comm in commodities])

set_odpairs(cmodel::AbstractConjugateSolver, ks::AbstractVector{Int}) =
    set_odpairs(cmodel, problem(cmodel).K[ks])

function set_demands(cmodel::AbstractConjugateSolver, demands)
    set_demands(cmodel, demands, tolled_arcs(problem(cmodel)))
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
