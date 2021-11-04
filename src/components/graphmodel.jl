## Build model
function make_graph_model(prob::Problem; silent=true)
    model = Model(() -> Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), silent)
    set_optimizer_attribute(model, MOI.NumberOfThreads(), 1)
    
    nv = prob.V
    na = length(prob.A)
    a1 = tolled_arcs(prob)
    model[:prob] = prob
    model[:a1] = a1

    c = cost_vector(prob)
    A = incidence_matrix(prob)

    @variable(model, 0 ≤ x[a=1:na])

    @objective(model, Min, sum(c .* x))

    @constraint(model, primalfeas, A * x .== zeros(nv))
    
    return model
end

# Set origin/destination
function set_OD!(model, orig, dest)
    prob = model[:prob]
    b = sourcesink_vector(prob, orig, dest)
    set_normalized_rhs.(model[:primalfeas], collect(b))
    return model
end

# Set fixed arcs
function set_fixedarcs!(model, fixed1, fixed0)
    a1 = model[:a1]
    x = model[:x]

    # Unfix everything
    for a in a1
        set_lower_bound(x[a], 0)
        set_upper_bound(x[a], Inf)
    end

    # Fix 1
    for a in fixed1
        set_lower_bound(x[a], 1)
        set_upper_bound(x[a], 1)
    end

    # Fix 0
    for a in fixed0
        set_upper_bound(x[a], 0)
    end

    return model
end