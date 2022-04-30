#=
A DualRepresentation must define:
- dualobj(dual): the follower's objective term of the dual (e.g. b'λ)
- formulate_dual!(model, dual): add the dual to the model
=#
abstract type DualRepresentation end

problem(dual::DualRepresentation) = dual.prob
dualobj(dual::DualRepresentation) = dual.dualobj

# Dual arc
mutable struct DualArc <: DualRepresentation
    prob::AbstractCommodityProblem
    λ::Union{Nothing,Vector{VariableRef}}
    dualobj::Union{Nothing,AffExpr}
    
    function DualArc(prob::AbstractCommodityProblem)
        return new(prob, nothing, nothing)
    end
end

function formulate_dual!(model::Model, dual::DualArc)
    prob = problem(dual)
    nv = nodes(prob)
    k = index(prob)
    
    c = cost_vector(prob)
    A = incidence_matrix(prob)
    b = sourcesink_vector(prob)
    
    dual.λ = λ = @variable(model, [1:nv], lower_bound = 0, base_name="λ[$k]")
    fix_var(λ[dest(prob)])
    
    t = remap_t(model, prob)
    tfull = expand_t(t, prob)
    
    @constraint(model, A' * λ .≤ c + tfull)
    dual.dualobj = b' * λ

    return
end

# Dual path
mutable struct DualPath <: DualRepresentation
    prob::PathPreprocessedProblem
    L::Union{Nothing,VariableRef}
    dualobj::Union{Nothing,AffExpr}

    function DualPath(prob::PathPreprocessedProblem)
        return new(prob, nothing, nothing)
    end
end

function formulate_dual!(model::Model, dual::DualPath)
    prob = problem(dual)
    k = index(prob)

    c = cost_vector(prob)
    δ = path_arc_incident_matrix(prob)

    dual.L = L = @variable(model, lower_bound = 0, base_name="L[$k]")

    t = remap_t(model, prob)
    tfull = expand_t(t, prob)

    @constraint(model, L .<= δ' * (c + tfull))
    dual.dualobj = L

    return
end
