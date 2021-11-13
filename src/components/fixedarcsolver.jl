function solve_fixedarcs_hungarian(graph::AbstractGraph, orig, dest, prob, fixed1, fixed0=[]; relaxed=false, tolls=zero_tolls(prob), setup=true)
    sources = [orig, [a.dst for a in prob.A[fixed1]]...]
    sinks = [dest, [a.src for a in prob.A[fixed1]]...]

    # Setting up
    if setup
        reset!(graph, prob)
        set_tolls!(graph, prob, tolls)
        disable_arcs!(graph, prob, fixed0)
        if !relaxed
            disable_arcs!(graph, prob, fixed1)
        end
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

    return cost + (isempty(fixed1) ? 0 : sum(prob.A[a].cost for a in fixed1))
end

function solve_fixedarcs_gurobi(model, orig, dest, fixed1, fixed0=[])
    set_OD!(model, orig, dest)
    set_fixedarcs!(model, fixed1, fixed0)
    optimize!(model)

    return termination_status(model) == MOI.OPTIMAL ? objective_value(model) : Inf
end

function _solve_fixedarcs_recur(graph::AbstractGraph, orig, dest, prob, fixed1, fixed0, depth, lastcost; tolls=zero_tolls(prob))
    if depth == 0 || isempty(fixed1)
        return solve_fixedarcs_hungarian(graph, orig, dest, prob, fixed1, fixed0, tolls=tolls, setup=false) + lastcost
    else
        mincost = Inf
        for a in fixed1
            arc = prob.A[a]
            _, cost = shortest_path(graph, orig, arc.src)
            
            thiscost = _solve_fixedarcs_recur(graph, arc.dst, dest, prob, setdiff(fixed1, a), fixed0, depth-1, lastcost + cost + arc.cost + tolls[a])
            if thiscost < mincost
                mincost = thiscost
            end
        end
        return mincost
    end
end

function solve_fixedarcs_recur(graph::AbstractGraph, orig, dest, prob, fixed1, fixed0; depth=2, tolls=zero_tolls(prob))
    reset!(graph, prob)
    set_tolls!(graph, prob, tolls)
    disable_arcs!(graph, prob, fixed0)
    disable_arcs!(graph, prob, fixed1)

    return _solve_fixedarcs_recur(graph, orig, dest, prob, fixed1, fixed0, depth, 0.0, tolls=tolls)
end