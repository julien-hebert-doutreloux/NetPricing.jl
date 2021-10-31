struct FeasiblePoint
    tolls::TollsDict                # Dict of tolls
    cost::Float64                   # Follower cost with set tolls
    rev::Float64                    # Leader revenue with set tolls
    paths::Vector{Vector{Int}}      # Optimal follower paths with set tolls
end

function Base.show(io::IO, point::FeasiblePoint)
    sorted_vals = round.(last.(sort(collect(point.tolls), by=first)), sigdigits=4)
    rev = round(point.rev, sigdigits=4)
    print(io, "Solution $sorted_vals with obj $rev")
end

## Get follower cost, revenue, and optiomal paths from tolls
function cost_revenue_paths(prob::Problem, tolls::AbstractTollsDict; graph = build_graph(prob), discount=1-eps())
    # Discount based on percentage to incentivize followers to use tolled arcs
    discounted_tolls = Dict(zip(keys(tolls), values(tolls) .* discount))
    set_tolls!(graph, prob, discounted_tolls)

    arctoindex = srcdst_to_index(prob)

    cost = 0.0
    rev = 0.0
    all_paths = Vector{Int}[]

    for comm in prob.K
        path = shortest_path(graph, comm.orig, comm.dest)[1]
        patharcs = tuple.(path[1:end-1], path[2:end])

        pathcost = 0.0
        pathrev = 0.0

        for arc in patharcs
            index = arctoindex[arc]
            pathcost += prob.A[index].cost
            if prob.A[index].toll
                pathrev += tolls[index]
            end
        end
        
        pathcost += pathrev
        
        cost += pathcost * comm.demand
        rev += pathrev * comm.demand
        push!(all_paths, path)
    end
    return cost, rev, all_paths
end

follower_cost(prob::Problem, tolls::AbstractTollsDict; kwargs...) = cost_revenue_paths(prob, tolls; kwargs...)[1]
revenue(prob::Problem, tolls::AbstractTollsDict; kwargs...) = cost_revenue_paths(prob, tolls; kwargs...)[2]
function base_cost(prob::Problem, tolls::AbstractTollsDict; kwargs...)
    cost, rev, _ = cost_revenue_paths(prob, tolls; kwargs...)
    return cost - rev
end

# Lower bound from set tolls
function lowerbound_from_tolls(prob::Problem, tolls::AbstractTollsDict; kwargs...)
    cost, rev, paths = cost_revenue_paths(prob, tolls; kwargs...)
    return FeasiblePoint(tolls, cost, rev, paths)
end

## Lower bound from solving the exact model
function lowerbound_from_exact_model(model)
    optimize!(model)
    prob = model[:prob]
    na = length(prob.A)
    nk = length(prob.K)

    rev = objective_value(model)
    # Prefer the dual obj because it's shorter
    cost = value(sum(model[:dualobj] .* demand_vector(prob)))

    # Extract tolls
    tvals = value_t(model)
    a1 = model[:a1]
    tolls = Dict(a => tvals[a] for a in a1)

    # Extract paths
    xvals = value_x(model)
    paths = Vector{Int}[]

    for k in 1:nk
        # Indices of all used arcs
        indices = [i for i in 1:na if xvals[k,i] > 0.5]
        # Build a src-dst dict
        srcdstdict = Dict(arc.src => arc.dst for arc in prob.A[indices])
        # Trace the path from orig to dest
        orig = prob.K[k].orig
        dest = prob.K[k].dest

        curr = orig
        path = Int[orig]
        while curr != dest
            curr = srcdstdict[curr]
            push!(path, curr)
        end

        push!(paths, path)
    end

    return FeasiblePoint(tolls, cost, rev, paths)
end