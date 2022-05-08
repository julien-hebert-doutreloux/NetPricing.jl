mutable struct ConjugateDynamicLinearModel <: AbstractConjugateSolver
    model::Union{Model,Nothing}
    prob::Problem
    num_commodities::Int
    weights::Vector{Float64}
    kwargs
    graph::AbstractGraph
    odpairs::Vector{Tuple{Int,Int}}
    tollfreecosts::Vector{Float64}
    nulltollcosts::Vector{Float64}
    nulltollarcs::Vector{BitSet}
    demands::Dict{Int,Float64}
end

# Interface
problem(cmodel::ConjugateDynamicLinearModel) = cmodel.prob
num_commodities(cmodel::ConjugateDynamicLinearModel) = cmodel.num_commodities
weights(cmodel::ConjugateDynamicLinearModel) = cmodel.weights

# Init a conjugate model to test bilevel feasibility between some paths, weights = 1
function ConjugateDynamicLinearModel(prob::Problem, num_commodities,
    weights=ones(num_commodities); kwargs...)

    graph = build_graph(prob)
    return ConjugateDynamicLinearModel(nothing, prob, num_commodities, weights, kwargs,
        graph, [], [], [], [], Dict())
end

# Init a conjugate model for the whole problem, weights = commodity demands, odpairs already set
function ConjugateDynamicLinearModel(prob::Problem; kwargs...)
    cmodel = ConjugateDynamicLinearModel(prob, length(prob.K), [comm.demand for comm in prob.K]; kwargs...)
    set_odpairs(cmodel, [(comm.orig, comm.dest) for comm in prob.K])
    return cmodel
end

# Set OD pairs
function set_odpairs(cmodel::ConjugateDynamicLinearModel, odpairs::AbstractVector{Tuple{Int,Int}})
    cmodel.odpairs = odpairs

    # Cache toll-free and null-toll paths
    @unpack prob, num_commodities, graph, tollfreecosts, nulltollcosts, nulltollarcs = cmodel
    empty!(tollfreecosts)
    empty!(nulltollcosts)
    empty!(nulltollarcs)

    nk = num_commodities
    a1set = BitSet(tolled_arcs(prob))
    arcdict = srcdst_to_index(prob)

    disable_tolls!(graph, prob)
    for k in 1:nk
        _, cost = shortest_path(graph, odpairs[k]...)
        push!(tollfreecosts, cost)
    end
    
    reset_tolls!(graph, prob)
    for k in 1:nk
        path, cost = shortest_path(graph, odpairs[k]...)
        path_a1 = BitSet(path_tolled_arcs(path, arcdict, a1set))
        push!(nulltollcosts, cost)
        push!(nulltollarcs, path_a1)
    end

    return nothing
end

# Set demands
function set_demands(cmodel::ConjugateDynamicLinearModel, demands, priority; discount=1-1e-5)
    a1 = tolled_arcs(cmodel.prob)
    cmodel.demands = Dict(a => get(demands, a, 0.) * (a ∈ priority ? discount : 1) for a in a1)
    return nothing
end

# Optimize
JuMP.objective_value(cmodel::ConjugateDynamicLinearModel) = objective_value(cmodel.model)
tvals(cmodel::ConjugateDynamicLinearModel) = value.(cmodel.model[:t]).data

function JuMP.optimize!(cmodel::ConjugateDynamicLinearModel)
    @unpack prob, graph, odpairs, nulltollarcs = cmodel
    nk = cmodel.num_commodities

    a1 = tolled_arcs(prob)
    a1set = BitSet(a1)
    arcdict = srcdst_to_index(prob)

    model = cmodel.model = _make_model(cmodel)
    L, t = model[:L], model[:t]

    # Start solving
    last_t = zeros(length(a1))
    last_path_a1 = copy(nulltollarcs)

    for i in 1:100
        optimize!(model)

        # Check if t is not changed
        tvals = value.(t).data
        all(isapprox.(last_t, tvals)) && break

        last_tvals = Dict(a1 .=> last_t)
        current_tvals = Dict(a1 .=> tvals)
        last_t = tvals

        # Set the current tolls to the graph
        set_tolls!(graph, prob, current_tvals)

        # Find the shortest paths
        for k in 1:nk
            # Skip empty commodities
            isempty(nulltollarcs[k]) && continue
            
            # Check if we can skip this commodity based on monotonicity:
            # For all a, a is used and t[a] decreases, or a is not used and t[a] increases
            all(a ->
                (a ∈ last_path_a1[k] && current_tvals[a] <= last_tvals[a]) ||
                (a ∉ last_path_a1[k] && current_tvals[a] >= last_tvals[a]),
                a1) && continue

            # If we can't skip, then find the shortest path and add a constraint
            path, cost = shortest_path(graph, odpairs[k]...)

            path_a1 = BitSet(path_tolled_arcs(path, arcdict, a1set))
            isempty(path_a1) && continue            # Toll-free path: already added
            last_path_a1[k] == path_a1 && continue  # Same as last path: no new constraint
            last_path_a1[k] = path_a1
            basecost = cost - sum(current_tvals[a] for a in path_a1)

            @constraint(model, L[k] ≤ basecost + sum(t[a] for a in path_a1))
        end

        # Fail-safe
        @assert(i < 100, "Too many iterations in optimize!(::ConjugateDynamicLinearModel)")
    end

    return
end

function _make_model(cmodel::ConjugateDynamicLinearModel)
    @unpack prob, num_commodities, weights, kwargs, graph, odpairs, tollfreecosts, nulltollcosts, nulltollarcs, demands = cmodel

    # Since we create a new model everytime, making a direct model is faster
    model = direct_model(Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), get(kwargs, :silent, true))
    set_optimizer_attribute(model, MOI.NumberOfThreads(), get(kwargs, :threads, nothing))

    nk = num_commodities
    a1 = tolled_arcs(prob)

    @variable(model, 0 ≤ L[k=1:nk] ≤ tollfreecosts[k])
    @variable(model, t[a=a1] ≥ 0)

    @objective(model, Max, weights' * L - sum(demands[a] * t[a] for a in a1))

    # Null-toll paths
    @constraint(model, [k=1:nk; !isempty(nulltollarcs[k])], L[k] ≤ nulltollcosts[k] + sum(t[a] for a in nulltollarcs[k]))

    return model
end
