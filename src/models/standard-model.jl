const _standard_model_refs = [
    :x,
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
    :dualobj
]

## Model for commodity
function add_standard_model!(model::Model, prob::AbstractCommodityProblem, M, N; sdtol=1e-10, linearize=true, dualanchor=true)
    nv = nodes(prob)
    na = length(arcs(prob))
    a1 = tolled_arcs(prob)
    k = index(prob)
    Amap = arcmap(prob)

    c = cost_vector(prob)
    A = incidence_matrix(prob)
    b = sourcesink_vector(prob)

    # Variables
    x = @variable(model, [a=1:na], lower_bound = 0, upper_bound = 1, binary = a in a1, base_name="x[$k]")
    λ = @variable(model, [i=1:nv], base_name="λ[$k]")
    tx = linearize ? @variable(model, [a=a1], lower_bound = 0, base_name="tx[$k]") : []
    t = JuMP.Containers.DenseAxisArray(model[:t][Amap[a1]].data, a1)

    sumtx = linearize ? sum(tx) : sum(t .* x[a1])

    a1dict = Dict(a => i for (i, a) in enumerate(tolled_arcs(parent(prob))))
    mapped_a1 = [a1dict[a] for a in Amap[a1]]
    M = @view M[mapped_a1]
    N = @view N[mapped_a1]

    # Expressions
    tfull = Array{Any}(zeros(na))
    for a in a1
        tfull[a] = t[a]
    end
    
    primalobj = c' * x + sumtx
    dualobj = b' * λ

    # Objective
    @objective(model, Max, objective_function(model) + demand(prob) * sumtx)

    # Constraints
    primalfeas = @constraint(model, A * x .== b)
    dualfeas = @constraint(model, A' * λ .≤ c + tfull)
    strongdual = @constraint(model, primalobj ≤ dualobj + sdtol)

    bilinear1 = linearize ? @constraint(model, tx .≥ 0) : []
    bilinear2 = linearize ? @constraint(model, tx .≤ M .* x[a1]) : []
    bilinear3 = linearize ? @constraint(model, t .- tx .≥ 0) : []
    bilinear4 = linearize ? @constraint(model, t .- tx .≤ N .* (1 .- x[a1])) : []

    # Dual anchor
    if dualanchor
        fix_var(λ[dest(prob)])
    end
    
    # References
    @pushrefs model _standard_model_refs

    return model
end

function add_standard_model!(model::Model, ::EmptyProblem, M, N; kwargs...)
    @pushemptyrefs model _standard_model_refs
    return model
end

## Build model
function standard_model(probs::Vector{<:AbstractCommodityProblem}; silent=false, threads=nothing, kwargs...)
    model = Model(() -> Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), silent)
    set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)

    # Big M
    parentprob = parent(probs[1])
    M, N = calculate_bigM(parentprob)

    a1 = tolled_arcs(parentprob)
    a1dict = Dict(a => i for (i, a) in enumerate(a1))

    @variable(model, 0 ≤ t[a=a1], upper_bound = N[a1dict[a]])
    @objective(model, Max, 0)

    # References
    @makerefs model _standard_model_refs

    # Add commodity
    for pprob in probs
        add_standard_model!(model, pprob, @view(M[index(pprob),:]), N; kwargs...)
    end
    
    return model
end

standard_model(prob::Problem; maxpaths=1000, kwargs...) = standard_model([preprocess(prob, k, maxpaths=maxpaths) for k in 1:length(prob.K)]; kwargs...)