function solve_fixedarcs_hungarian(graph::AbstractGraph, orig, dest, prob, fixed1, fixed0=[]; relaxed=false)
    sources = [orig, [a.dst for a in prob.A[fixed1]]...]
    sinks = [dest, [a.src for a in prob.A[fixed1]]...]

    # Setting up
    reset!(graph, prob)
    disable_arcs!(graph, prob, fixed0)

    if !relaxed
        disable_arcs!(graph, prob, fixed1)
    end

    # Filling the cost matrix
    distances = [dijkstra_shortest_paths(graph, src).dists for src in sources]
    costs = [distances[row][sink] for row in 1:length(distances), sink in sinks]

    # If any row or column is Inf, return Inf
    # Hungarian.jl stucks if this is not checked
    if any(isinf, minimum(costs, dims=1)) || any(isinf, minimum(costs, dims=2))
        return Inf
    end

    # Solve the assignment
    assignment, cost = hungarian(costs)

    return cost + sum(prob.A[a].cost for a in fixed1)
end

function solve_fixedarcs_gurobi(model, orig, dest, fixed1, fixed0=[])
    set_OD!(model, orig, dest)
    set_fixedarcs!(model, fixed1, fixed0)
    optimize!(model)

    return termination_status(model) == MOI.OPTIMAL ? objective_value(model) : Inf
end