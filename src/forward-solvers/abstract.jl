abstract type AbstractForwardSolver end

# Interface
function set_tolls(fmodel::AbstractForwardSolver, tolls)
    set_tolls(fmodel, tolls, tolled_arcs(problem(fmodel)))
end

# Solve for f(t) and w
function solve(fmodel::AbstractForwardSolver)
    optimize!(fmodel)
    return objective_value(fmodel), wvals(fmodel)
end

function solve(fmodel::AbstractForwardSolver, tolls, args...)
    set_tolls(fmodel, tolls, args...)
    return solve(fmodel)
end
