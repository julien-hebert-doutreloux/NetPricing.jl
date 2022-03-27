# Dual arc
function formulate_dual(model::Model, form::GeneralFormulation{<:Any,DualArc})
    prob = problem(dual(form))
    nv = nodes(prob)
    k = index(prob)
    
    c = cost_vector(prob)
    A = incidence_matrix(prob)
    b = sourcesink_vector(prob)
    
    λ = @variable(model, [1:nv], lower_bound = 0, base_name="λ[$k]")
    fix_var(λ[dest(prob)])
    
    t = remap_t(model, prob)
    tfull = expand_t(t, prob)
    
    @constraint(model, A' * λ .≤ c + tfull)
    dualobj = b' * λ

    return dualobj
end

# Dual path
function formulate_dual(model::Model, form::GeneralFormulation{<:Any,DualPath})
    prob = problem(dual(form))
    k = index(prob)

    c = cost_vector(prob)
    δ = path_arc_incident_matrix(prob)

    L = @variable(model, lower_bound = 0, base_name="L[$k]")

    t = remap_t(model, prob)
    tfull = expand_t(t, prob)

    @constraint(model, L .<= δ' * (c + tfull))
    dualobj = L

    return dualobj
end
