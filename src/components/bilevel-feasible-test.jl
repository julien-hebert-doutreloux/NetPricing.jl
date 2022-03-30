function bilevel_feasible_test_model(prob::Problem, num_commodities::Int; silent=true, threads=nothing)
    model = Model(() -> Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), silent)
    set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)

    nk = num_commodities
    na = length(arcs(prob))
    a1 = tolled_arcs(prob)

    c = cost_vector(prob)
    A = incidence_matrix(prob)

    @variable(model, 0 ≤ x[a=1:na, k=1:nk])
    @objective(model, Min, sum(c' * x))
    @constraint(model, networkflow, A * x .== 0)
    @constraint(model, demandlimit[a=a1], sum(x[a,k] for k in 1:nk) ≤ nk)

    model[:prob] = prob
    model[:nk] = nk

    return model
end

function set_bilevel_feasible_test_odpairs(model::Model, odpairs::AbstractVector{Tuple{Int,Int}})
    prob = model[:prob]
    b = hcat([sourcesink_vector(prob, o, d) for (o, d) in odpairs]...)
    set_normalized_rhs.(model[:networkflow], b)
    return model
end

set_bilevel_feasible_test_odpairs(model::Model, commodities::AbstractVector{Commodity}) =
    set_bilevel_feasible_test_odpairs(model, [(comm.orig, comm.dest) for comm in commodities])

set_bilevel_feasible_test_odpairs(model::Model, ks::AbstractVector{Int}) =
    set_bilevel_feasible_test_odpairs(model, model[:prob].K[ks])

function is_bilevel_feasible(model::Model, paths; set_odpairs=true)
    prob = model[:prob]
    a1 = tolled_arcs(prob)

    # Set od pairs if required
    if set_odpairs
        set_bilevel_feasible_test_odpairs(model, [(path[1], path[end]) for path in paths])
    end

    # Extract the set of tolled arcs for each path
    arcdict = srcdst_to_index(prob)
    a1set = BitSet(a1)
    tolled_sets = BitSet[]
    for path in paths
        tolled_set = BitSet(path_tolled_arcs(path, arcdict, a1set))
        push!(tolled_sets, tolled_set)
    end

    # Set the demand
    for a in a1
        demand = count(set -> a in set, tolled_sets)
        set_normalized_rhs(model[:demandlimit][a], demand)
    end

    # Solve the model
    optimize!(model)

    # Check the objective value
    objval = objective_value(model)
    arccosts = srcdst_to_cost(prob)
    expecting = sum(get_path_cost(path, arccosts) for path in paths)

    return isapprox(objval, expecting, rtol = 1e-6)
end
is_bilevel_feasible(prob::Problem, paths) = is_bilevel_feasible(bilevel_feasible_test_model(prob, length(paths)), paths)
