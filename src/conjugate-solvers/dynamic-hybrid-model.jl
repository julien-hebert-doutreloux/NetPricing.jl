# Similar to ConjugateDynamicHybridModel but using ForwardHybridSolver
mutable struct ConjugateDynamicHybridModel <: AbstractConjugateSolver
    model::Union{Model,Nothing}
    probs::Vector{AbstractCommodityProblem}
    weights::Vector{Float64}
    kwargs
    fsolver::ForwardHybridSolver
    demands::Dict{Int,Float64}
end

# Interface
problem(cmodel::ConjugateDynamicHybridModel) = cmodel.probs |> first |> parent
num_commodities(cmodel::ConjugateDynamicHybridModel) = length(problem(cmodel).K)
weights(cmodel::ConjugateDynamicHybridModel) = cmodel.weights

# Init a conjugate model for the whole problem, weights = commodity demands, odpairs already set
function ConjugateDynamicHybridModel(probs::Vector{<:AbstractCommodityProblem}, weights=demand.(probs); maxpaths=1000, kwargs...)
    fsolver = ForwardHybridSolver(probs; maxpaths)
    return ConjugateDynamicHybridModel(nothing, probs, weights, kwargs, fsolver, Dict())
end

# Set demands
function set_demands(cmodel::ConjugateDynamicHybridModel, demands, priority; discount=default_discount())
    a1 = cmodel |> problem |> tolled_arcs
    cmodel.demands = Dict(a => get(demands, a, 0.) * (a ∈ priority ? discount : 1) for a in a1)
    return nothing
end

# Optimize
JuMP.objective_value(cmodel::ConjugateDynamicHybridModel) = objective_value(cmodel.model)
tvals(cmodel::ConjugateDynamicHybridModel) = value.(cmodel.model[:t]).data

function JuMP.optimize!(cmodel::ConjugateDynamicHybridModel)
    @unpack fsolver = cmodel

    nk = num_commodities(cmodel)
    a1 = tolled_arcs(problem(cmodel))

    model = cmodel.model = _make_model(cmodel)
    L, t = model[:L], model[:t]

    # Start solving
    last_tvals = Dict(a => Inf for a in a1)
    last_path_a1 = [BitSet() for _ in 1:nk]

    for i in 1:10000
        optimize!(model)

        # Check if t is not changed
        tvals = value.(t).data
        all(isapprox(last_tvals[a], t) for (a, t) in zip(a1, tvals)) && break
        current_tvals = Dict(zip(a1, tvals))

        # Find the shortest paths
        for (k, solver) in enumerate(fsolver.solvers)
            # Skip empty commodities
            if solver isa CommodityForwardEmptySolver
                continue
            end

            # Check if we can skip this commodity based on monotonicity:
            # For all a, a is used and t[a] decreases, or a is not used and t[a] increases
            all(a ->
                (a ∈ last_path_a1[k] && current_tvals[a] <= last_tvals[a]) ||
                (a ∉ last_path_a1[k] && current_tvals[a] >= last_tvals[a]),
                a1) && continue

            # If we can't skip, solve the forward model
            set_tolls(solver, current_tvals)
            optimize!(solver)
            cost = objective_value(solver)
            path_a1 = optimal_tolled_set(solver)

            isempty(path_a1) && continue            # Toll-free path: already added
            last_path_a1[k] == path_a1 && continue  # Same as last path: no new constraint
            last_path_a1[k] = path_a1
            basecost = cost - sum(current_tvals[a] for a in path_a1)

            @constraint(model, L[k] ≤ basecost + sum(t[a] for a in path_a1))
        end

        last_tvals = current_tvals

        # Fail-safe
        @assert(i < 10000, "Too many iterations in optimize!(::ConjugateDynamicHybridModel)")
    end

    return
end

function _make_model(cmodel::ConjugateDynamicHybridModel)
    @unpack weights, kwargs, fsolver, demands = cmodel

    # Since we create a new model everytime, making a direct model is faster
    model = direct_model(Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), get(kwargs, :silent, true))
    set_optimizer_attribute(model, MOI.NumberOfThreads(), get(kwargs, :threads, nothing))

    nk = num_commodities(cmodel)
    a1 = tolled_arcs(problem(cmodel))

    solve(fsolver, Dict(a1 .=> Inf))
    tollfreecosts = objective_value.(fsolver.solvers)

    @variable(model, 0 ≤ L[k=1:nk] ≤ tollfreecosts[k])
    @variable(model, t[a=a1] ≥ 0)

    @objective(model, Max, weights' * L - sum(demands[a] * t[a] for a in a1))

    return model
end
