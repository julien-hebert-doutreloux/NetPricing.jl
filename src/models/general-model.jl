const _general_model_refs = [
    :x,
    :z,
    :λ,
    :tx,
    :primalfeas,
    :dualfeas,
    :strongdual,
    :bilinear1,
    :bilinear2,
    :bilinear3,
    :bilinear4,
    :primalobj,
    :dualobj,
    :conversion
]

## Model for commodity
function add_general_model!(model::Model, ::EmptyProblem, args...; kwargs...)
    @pushemptyrefs model _general_model_refs
    return model
end

function add_general_model!(model::Model, prob::AbstractCommodityProblem, primal_repr, dual_repr, M, N;
    sdtol=1e-10,            # Strong duality tolerance
    dualanchor=true,        # Fix λ_d to 0
    dualbound=true,         # Add lower bound to λ
    )

    x, primalfeas, primalobj, z, conversion = primal_repr(model, prob)
    λ, dualfeas, dualobj = dual_repr(model, prob; dualanchor=dualanchor, dualbound=dualbound)
    tx, sumtx, bilinear1, bilinear2, bilinear3, bilinear4 = linearization(model, prob, x, M, N)

    strongdual = @constraint(model, primalobj + sumtx ≤ dualobj + sdtol)

    # Objective
    @objective(model, Max, objective_function(model) + demand(prob) * sumtx)
    
    # References
    @pushrefs model _general_model_refs

    return model
end

## Build model
function general_model(probs::Vector{<:AbstractCommodityProblem}, primal_repr, dual_repr; silent=false, threads=nothing, bigM_maxpaths=100, kwargs...)
    model = Model(() -> Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), silent)
    set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)

    # Big M
    parentprob = parent(probs[1])
    _, M, N = calculate_bigM_paths(probs, threads=threads, maxpaths=bigM_maxpaths)

    a1 = tolled_arcs(parentprob)
    a1dict = Dict(a => i for (i, a) in enumerate(a1))

    @variable(model, 0 ≤ t[a=a1], upper_bound = N[a1dict[a]])
    @objective(model, Max, 0)

    # References
    @makerefs model _general_model_refs

    # Add commodity
    for pprob in probs
        add_general_model!(model, pprob, primal_repr, dual_repr, @view(M[index(pprob),:]), N; kwargs...)
    end

    model[:probs] = probs
    
    return model
end

general_model(prob::Problem, primal_repr, dual_repr; maxpaths=1000, kwargs...) =
    general_model([preprocess(prob, k, maxpaths=maxpaths) for k in 1:length(prob.K)], primal_repr, dual_repr; kwargs...)