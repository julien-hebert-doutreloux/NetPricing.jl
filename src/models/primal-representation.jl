#=
A PrimalRepresentation must define:
- primalobj(primal): the follower's objective term of the primal (e.g. c'x)
- formulate_primal!(model, primal): add the primal to the model
=#
abstract type PrimalRepresentation end

problem(primal::PrimalRepresentation) = primal.prob
primalobj(primal::PrimalRepresentation) = primal.primalobj

# Primal arc
mutable struct PrimalArc <: PrimalRepresentation
    prob::AbstractCommodityProblem
    x::Union{Nothing,Vector{VariableRef}}
    primalobj::Union{Nothing,AffExpr}

    function PrimalArc(prob::AbstractCommodityProblem)
        return new(prob, nothing, nothing)
    end
end

function formulate_primal!(model::Model, primal::PrimalArc)
    prob = problem(primal)
    na = length(arcs(prob))
    k = index(prob)
    
    c = cost_vector(prob)
    A = incidence_matrix(prob)
    b = sourcesink_vector(prob)

    primal.x = x = @variable(model, [1:na], lower_bound = 0, upper_bound = 1, base_name="x[$k]")

    @constraint(model, A * x .== b)
    primal.primalobj = c' * x

    return
end

# Primal path
mutable struct PrimalPath <: PrimalRepresentation
    prob::PathPreprocessedProblem
    z::Union{Nothing,Vector{VariableRef}}
    x::Union{Nothing,DenseAxisArray{VariableRef}}
    primalobj::Union{Nothing,AffExpr}

    function PrimalPath(prob::PathPreprocessedProblem)
        return new(prob, nothing, nothing, nothing)
    end
end


function formulate_primal!(model::Model, primal::PrimalPath)
    prob = problem(primal)
    all_paths = prob.paths
    np = length(all_paths)
    a1 = tolled_arcs(prob)
    k = index(prob)

    c = cost_vector(prob)
    δ = path_arc_incident_matrix(prob)
    δ1 = @view δ[a1, :]

    primal.z = z = @variable(model, [1:np], lower_bound = 0, upper_bound = 1, base_name="z[$k]")
    primal.x = x = @variable(model, [a1], lower_bound = 0, upper_bound = 1, base_name="x[$k]")

    @constraint(model, sum(z) == 1)
    @constraint(model, x .== δ1 * z)
    primal.primalobj = c' * δ * z

    return
end
