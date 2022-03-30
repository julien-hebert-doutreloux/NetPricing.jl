# Representations
abstract type DualRepresentation end

problem(dual::DualRepresentation) = dual.prob

mutable struct DualArc <: DualRepresentation
    prob::AbstractCommodityProblem
    λ::Union{Nothing,Vector{VariableRef}}
    
    function DualArc(prob::AbstractCommodityProblem)
        return new(prob, nothing)
    end
end

mutable struct DualPath <: DualRepresentation
    prob::PathPreprocessedProblem
    L::Union{Nothing,VariableRef}

    function DualPath(prob::PathPreprocessedProblem)
        return new(prob, nothing)
    end
end

# Dual arc
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
    dualobj = b' * λ

    return dualobj
end

# Dual path
function formulate_dual!(model::Model, dual::DualPath)
    prob = problem(dual)
    k = index(prob)

    c = cost_vector(prob)
    δ = path_arc_incident_matrix(prob)

    dual.L = L = @variable(model, lower_bound = 0, base_name="L[$k]")

    t = remap_t(model, prob)
    tfull = expand_t(t, prob)

    @constraint(model, L .<= δ' * (c + tfull))
    dualobj = L

    return dualobj
end
