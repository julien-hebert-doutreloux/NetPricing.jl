macro makerefs(model, args...)
    model = esc(model)
    ex = Expr(:block)
    ex.args = [:($model[$(Meta.quot(arg))] = []) for arg in args]
    return ex
end

macro pushrefs(model, args...)
    model = esc(model)
    ex = Expr(:block)
    ex.args = [:(push!($model[$(Meta.quot(arg))], $(esc(arg)))) for arg in args]
    return ex
end

macro pushemptyrefs(model, args...)
    model = esc(model)
    ex = Expr(:block)
    ex.args = [:(push!($model[$(Meta.quot(arg))], [])) for arg in args]
    return ex
end

## Model for commodity
function add_standard_model!(model::Model, prob::AbstractCommodityProblem, M, N; sdtol=1e-10)
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
    tx = @variable(model, [a=a1], lower_bound = 0, base_name="tx[$k]")
    t = JuMP.Containers.DenseAxisArray(model[:t][Amap[a1]].data, a1)

    a1dict = Dict(a => i for (i, a) in enumerate(tolled_arcs(parent(prob))))
    mapped_a1 = [a1dict[a] for a in Amap[a1]]
    M = @view M[mapped_a1]
    N = @view N[mapped_a1]

    # Expressions
    tfull = Array{Any}(zeros(na))
    for a in a1
        tfull[a] = t[a]
    end
    
    primalobj = c' * x + sum(tx)
    dualobj = b' * λ

    # Objective
    set_objective_coefficient.(model, tx, demand(prob))

    # Constraints
    primalfeas = @constraint(model, A * x .== b)
    dualfeas = @constraint(model, A' * λ .≤ c + tfull)
    strongdual = @constraint(model, primalobj ≤ dualobj + sdtol)
    bilinear1 = @constraint(model, tx .≥ 0)
    bilinear2 = @constraint(model, tx .≤ M .* x[a1])
    bilinear3 = @constraint(model, t .- tx .≥ 0)
    bilinear4 = @constraint(model, t .- tx .≤ N .* (1 .- x[a1]))
    
    # References
    @pushrefs model x λ tx primalfeas dualfeas strongdual bilinear1 bilinear2 bilinear3 bilinear4 primalobj dualobj

    return model
end

function add_standard_model!(model::Model, ::EmptyProblem, M, N; sdtol=1e-10)
    @pushemptyrefs model x λ tx primalfeas dualfeas strongdual bilinear1 bilinear2 bilinear3 bilinear4 primalobj dualobj
    return model
end

## Build model
function standard_model(prob::Problem; sdtol=1e-10, silent=false, threads=nothing, maxpaths=1000)
    model = Model(() -> Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), silent)
    set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)

    # Preprocess
    pprobs = preprocess(prob, maxpaths=maxpaths)

    # Big M
    M, N = calculate_bigM(prob)

    a1 = tolled_arcs(prob)
    @variable(model, t[a=a1] ≥ 0)
    @objective(model, Max, 0)

    # References
    @makerefs model x λ tx primalfeas dualfeas strongdual bilinear1 bilinear2 bilinear3 bilinear4 primalobj dualobj

    # Add commodity
    for pprob in pprobs
        add_standard_model!(model, pprob, @view(M[index(pprob),:]), N, sdtol=sdtol)
    end
    
    return model
end
