function inverse_model(prob::Problem, num_commodities::Int; silent=false, threads=nothing)
    model = Model(() -> Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), silent)
    set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)

    nk = num_commodities
    na = length(arcs(prob))
    nv = nodes(prob)
    a1 = tolled_arcs(prob)

    c = cost_vector(prob)
    A = incidence_matrix(prob)

    @variable(model, λ[i=1:nv, k=1:nk])
    @variable(model, t[a=a1] ≥ 0)

    @objective(model, Max, 0)

    tfull = Array{Any}(zeros(na))
    for a in a1
        tfull[a] = t[a]
    end
    @constraint(model, dualfeas, A' * λ .- tfull .≤ c)

    model[:prob] = prob
    model[:nk] = nk
    
    return model
end

# Set OD pairs
function set_inverse_model_odpairs(model::Model, odpairs::AbstractVector{Tuple{Int,Int}})
    prob = model[:prob]
    b = hcat([sourcesink_vector(prob, o, d) for (o, d) in odpairs]...)
    set_objective_coefficient.(model, model[:λ], b)
    return model
end

set_inverse_model_odpairs(model::Model, commodities::AbstractVector{Commodity}) =
    set_inverse_model_odpairs(model, [(comm.orig, comm.dest) for comm in commodities])

set_inverse_model_odpairs(model::Model, ks::AbstractVector{Int}) =
    set_inverse_model_odpairs(model, model[:prob].K[ks])

# Set paths
function set_inverse_model_paths(model::Model, paths, discount::AbstractDict{Int,Float64}; set_odpairs=true)
    prob = model[:prob]
    nk = model[:nk]
    a1 = tolled_arcs(prob)

    # Set od pairs if required
    if set_odpairs
        set_inverse_model_odpairs(model, [(path[1], path[end]) for path in paths])
    end

    # Extract the set of tolled arcs for each path
    arcdict = srcdst_to_index(prob)
    a1set = BitSet(a1)
    tolled_sets = [BitSet(path_tolled_arcs(path, arcdict, a1set)) for path in paths]

    # Set the demand
    for a in a1
        demand = count(set -> a in set, tolled_sets)
        set_objective_coefficient(model, model[:t][a], -demand * get(discount, a, 1))
    end

    return model
end

function set_inverse_model_paths(model::Model, paths, discount::Real = 1-1e-5; kwargs...)
    a1 = tolled_arcs(model[:prob])
    set_inverse_model_paths(model, paths, Dict(a1 .=> discount); kwargs...)
end
