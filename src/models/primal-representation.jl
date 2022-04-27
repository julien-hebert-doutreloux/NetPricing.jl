# Representations
abstract type PrimalRepresentation end

problem(primal::PrimalRepresentation) = primal.prob

mutable struct PrimalArc <: PrimalRepresentation
    prob::AbstractCommodityProblem
    x::Union{Nothing,Vector{VariableRef}}

    function PrimalArc(prob::AbstractCommodityProblem; kwargs...)
        return new(prob, nothing)
    end
end

mutable struct PrimalPath <: PrimalRepresentation
    prob::PathPreprocessedProblem
    binary_x::Bool
    z::Union{Nothing,Vector{VariableRef}}
    x::Union{Nothing,DenseAxisArray{VariableRef}}

    function PrimalPath(prob::PathPreprocessedProblem; binary_x::Bool=false)
        return new(prob, binary_x, nothing, nothing)
    end
end

# Primal arc
function formulate_primal!(model::Model, primal::PrimalArc)
    prob = problem(primal)
    na = length(arcs(prob))
    a1 = tolled_arcs(prob)
    k = index(prob)
    
    c = cost_vector(prob)
    A = incidence_matrix(prob)
    b = sourcesink_vector(prob)

    primal.x = x = @variable(model, [a=1:na], lower_bound = 0, upper_bound = 1, binary = a in a1, base_name="x[$k]")

    @constraint(model, A * x .== b)
    primalobj = c' * x

    return x, primalobj
end

# Primal path
function formulate_primal!(model::Model, primal::PrimalPath)
    prob = problem(primal)
    all_paths = prob.paths
    np = length(all_paths)
    a1 = tolled_arcs(prob)
    k = index(prob)

    c = cost_vector(prob)
    δ = path_arc_incident_matrix(prob)
    δ1 = @view δ[a1, :]

    primal.z = z = @variable(model, [1:np], lower_bound = 0, upper_bound = 1, binary = !primal.binary_x, base_name="z[$k]")
    primal.x = x = @variable(model, [a1], lower_bound = 0, upper_bound = 1, binary = primal.binary_x, base_name="x[$k]")

    @constraint(model, sum(z) == 1)
    @constraint(model, x .== δ1 * z)

    primalobj = c' * δ * z

    return x, primalobj
end

# Linearization
function linearization(model::Model, prob::AbstractCommodityProblem, x, M, N)
    parentprob = parent(prob)
    a1 = tolled_arcs(prob)
    k = index(prob)
    Amap = arcmap(prob)

    tx = @variable(model, [a=a1], lower_bound = 0, base_name="tx[$k]")
    t = remap_t(model, prob)

    sumtx = sum(tx)

    a1dict = Dict(a => i for (i, a) in enumerate(tolled_arcs(parentprob)))
    mapped_a1 = [a1dict[a] for a in Amap[a1]]
    M = @view M[mapped_a1]
    N = @view N[mapped_a1]

    @constraint(model, tx .≤ M .* x[a1])
    @constraint(model, t .- tx .≥ 0)
    @constraint(model, t .- tx .≤ N .* (1 .- x[a1]))

    return sumtx
end