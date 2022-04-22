struct ConjugateLinearModel <: AbstractConjugateSolver
    model::Model
    prob::Problem
    num_commodities::Int
    weights::Vector{Float64}
end

# Interface
problem(cmodel::ConjugateLinearModel) = cmodel.prob
num_commodities(cmodel::ConjugateLinearModel) = cmodel.num_commodities
weights(cmodel::ConjugateLinearModel) = cmodel.weights

# Init a conjugate model to test bilevel feasibility between some paths, weights = 1
function ConjugateLinearModel(prob::Problem, num_commodities,
    weights=ones(num_commodities);
    silent=true, threads=nothing)

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
    
    return ConjugateLinearModel(model, prob, num_commodities, weights)
end

# Init a conjugate model for the whole problem, weights = commodity demands, odpairs already set
function ConjugateLinearModel(prob::Problem; kwargs...)
    cmodel = ConjugateLinearModel(prob, length(prob.K), [comm.demand for comm in prob.K]; kwargs...)
    set_odpairs(cmodel, [(comm.orig, comm.dest) for comm in prob.K])
    return cmodel
end

# Set OD pairs
function set_odpairs(cmodel::ConjugateLinearModel, odpairs::AbstractVector{Tuple{Int,Int}})
    model, prob, η = cmodel.model, cmodel.prob, cmodel.weights
    b = hcat([sourcesink_vector(prob, o, d) for (o, d) in odpairs]...) .* η'
    set_objective_coefficient.(model, model[:λ], b)
    return nothing
end

# Set demands
function set_demands(cmodel::ConjugateLinearModel, demands, discount::AbstractDict{Int,<:Real})
    model = cmodel.model
    for (a, w) in demands
        set_objective_coefficient(model, model[:t][a], -w * get(discount, a, 1))
    end
    return nothing
end

# Optimize
JuMP.optimize!(cmodel::ConjugateLinearModel) = optimize!(cmodel.model)
JuMP.objective_value(cmodel::ConjugateLinearModel) = objective_value(cmodel.model)
tvals(cmodel::ConjugateLinearModel) = value.(cmodel.model[:t]).data
