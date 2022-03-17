# Dual arc
function dual_arc(model::Model, prob::AbstractCommodityProblem;
    dualanchor=true,        # Fix λ_d to 0
    dualbound=true,         # Add lower bound to λ
)
    nv = nodes(prob)
    na = length(arcs(prob))
    a1 = tolled_arcs(prob)
    k = index(prob)
    Amap = arcmap(prob)
    
    c = cost_vector(prob)
    A = incidence_matrix(prob)
    b = sourcesink_vector(prob)
    
    # Variables
    λ = @variable(model, [i=1:nv], base_name="λ[$k]")
    t = JuMP.Containers.DenseAxisArray(model[:t][Amap[a1]].data, a1)
    
    # Expressions
    tfull = Array{Any}(zeros(na))
    for a in a1
        tfull[a] = t[a]
    end
    
    # Constraints
    dualfeas = @constraint(model, A' * λ .≤ c + tfull)
    dualobj = b' * λ

    # Dual bounds
    if dualbound
        set_lower_bound.(λ, 0)
    end
    
    # Dual anchor
    if dualanchor
        fix_var(λ[dest(prob)])
    end

    return λ, dualfeas, dualobj
end
