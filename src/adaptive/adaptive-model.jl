## Build model
function adaptive_model(prob::Problem;
    silent=false,
    threads=nothing,
    timelimit=nothing,
    sdtol=1e-10,
    nonadaptive=false,
    dualanchor=true,
    preprocess=nothing)

    # model = Model(() -> Gurobi.Optimizer(current_env()))
    model = direct_model(Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), silent)
    set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)
    set_optimizer_attribute(model, MOI.TimeLimitSec(), timelimit)
    if !nonadaptive
        set_optimizer_attribute(model, "LazyConstraints", 1)
        set_optimizer_attribute(model, "Heuristics", 0)         # Turn off to avoid unnecessary graph expansion
    end

    # Parameters
    nv = prob.V
    na = length(prob.A)
    nk = length(prob.K)
    a1 = tolled_arcs(prob)

    model[:prob] = prob
    model[:a1] = a1
    model[:sdtol] = sdtol

    c = cost_vector(prob)
    A = incidence_matrix(prob)
    b = [sourcesink_vector(prob, k) for k in 1:nk]
    d = demand_vector(prob)

    Aout = max.(A, 0)
    Ain = max.(-A, 0)

    # Big M
    M, N = calculate_bigM(prob)
    odmax = calculate_odmax(prob)
    mincost = calculate_mincost_fixedarc(prob) * (1 - 1e-1)

    # Variables
    @variables(model,
        begin
            obj ≥ 0                                             # Leader objective
            0 ≤ x[k=1:nk, a=1:na] ≤ 1, (binary = a in a1)       # All arcs
            λ[k=1:nk, i=1:nv]                                   # Duals
            t[a=a1] ≥ 0                                         # Tolls
            tx[k=1:nk, a=a1] ≥ 0                                # Bilinear
            r[k=1:nk] ≥ 0                                       # Free profit
            τ[k=1:nk] ≥ 0                                       # Profit per commodity
        end)

    # Expressions
    tfull = Array{Any}(zeros(na))
    for a in a1
        tfull[a] = t[a]
    end
    
    model[:cx] = [c .* x[k,:] for k in 1:nk]
    dualobj = model[:dualobj] = [b[k]' * λ[k,:] for k in 1:nk]

    # Objective
    @objective(model, Max, obj)
    @constraints(model,
        begin
            objdef, obj == sum(d .* τ)
            taudef[k=1:nk], τ[k] == sum(tx[k,:]) + r[k]
        end)

    # Constraints
    @constraints(model,
        begin
            flowlimitout[k=1:nk], Aout * x[k,:] .≤ 1
            flowlimitin[k=1:nk], Ain * x[k,:] .≤ 1
            dualfeasdefault[k=1:nk], b[k]' * λ[k,:] ≤ odmax[k]
            strongdual[k=1:nk], c' * x[k,:] + τ[k] == dualobj[k] + sdtol
            bilinear1[k=1:nk], tx[k,:] .≥ 0
            bilinear2[k=1:nk], tx[k,:] .≤ M[k,:] .* x[k,a1]
            bilinear3[k=1:nk], t .- tx[k,:] .≥ 0
            bilinear4[k=1:nk], t .- tx[k,:] .≤ N .* (1 .- x[k,a1])
        end)

    # Anchor the dual at the destination
    if dualanchor
        for (k, comm) in enumerate(prob.K)
            fix_var(λ[k, comm.dest])
        end
    end

    # Build primalfeas and dualfeas beforehand (add to the model if non-adaptive)
    model[:primalfeas] = [@build_constraint(A * x[k,:] .== b[k]) for k in 1:nk]
    model[:dualfeas] = [@build_constraint(A' * λ[k,:] .≤ c + tfull) for k in 1:nk]
    model[:strongdualguide] = [@build_constraint(mincost[k,:] .* x[k,:] .+ τ[k] .≤ dualobj[k] .+ sdtol) for k in 1:nk]
    if nonadaptive
        for k in 1:nk
            add_constraint.(Ref(model), model[:primalfeas][k])
            add_constraint.(Ref(model), model[:dualfeas][k])
            add_constraint.(Ref(model), model[:strongdualguide][k])
        end
    end

    # Preprocessing
    isnothing(preprocess) && (preprocess = [UnprocessedProblem(prob, k) for k in 1:nk])
    model[:relevants] = [BitSet(used_nodes(pprob)) for pprob in preprocess]
    for k in 1:nk
        irrelevants_arcs = setdiff!(BitSet(1:na), used_arcs(preprocess[k]))
        fix_var.(x[k, collect(irrelevants_arcs)])
    end

    # Lazy constraint
    if !nonadaptive
        MOI.set(model, MOI.LazyConstraintCallback(), adaptive_lazy_routine(model; sdtol=sdtol))
    end
    
    return model
end

# Lazy constraint routine
function adaptive_lazy_routine(model; sdtol=1e-10)
    prob = model[:prob]

    A = incidence_matrix(prob)
    b = [sourcesink_vector(prob, k) for k in 1:length(prob.K)]

    a1 = tolled_arcs(prob)
    arcdict = srcdst_to_index(prob)

    relevants = model[:relevants]

    explicitnodes = model[:explicitnodes] = [BitSet() for _ in prob.K]
    explicitarcs = model[:explicitarcs] = [BitSet() for _ in prob.K]

    primalfeas = model[:primalfeas]
    dualfeas = model[:dualfeas]
    strongdualguide = model[:strongdualguide]

    t = model[:t]
    x = model[:x]

    dualgraph = build_graph(prob)

    return function adaptive_lazy_routine_inner(cb_data)
        status = callback_node_status(cb_data, model)
        (status == MOI.CALLBACK_NODE_STATUS_INTEGER) || return

        tvals = callback_value.(cb_data, t)
        set_tolls!(dualgraph, prob, Dict(a1 .=> tvals))

        primaloptimals = trues(length(prob.K))
        dualoptimals = trues(length(prob.K))
        
        # Add cuts
        for (k, comm) in enumerate(prob.K)
            orig, dest = comm.orig, comm.dest

            # Primal test: test flow unbalanced (only input > output)
            xvals = callback_value.(cb_data, x[k,:])
            unbalanced = findall(A * round.(Int, xvals) .!= b[k])
            filter!(i -> i ∈ relevants[k], unbalanced)

            if !isempty(unbalanced)
                # If not explicit, make it explicit
                primaloptimals[k] = false
                union!(explicitnodes[k], unbalanced)
                MOI.submit.(model, MOI.LazyConstraint(cb_data), primalfeas[k][unbalanced])

                # Add the guides
                # guide_arcs = findall(arc -> (arc.src ∈ explicitnodes[k]) ⊻ (arc.dst ∈ explicitnodes[k]), prob.A)
                # MOI.submit.(model, MOI.LazyConstraint(cb_data), strongdualguide[k][guide_arcs])
            end

            # Dual test: check the optimal path in the graph
            dualpath, _ = shortest_path(dualgraph, orig, dest)
            implicitarcs = Int[]
            for i in 1:(length(dualpath)-1)
                push!(implicitarcs, arcdict[(dualpath[i], dualpath[i+1])])
            end
            filter!(a -> a ∉ explicitarcs[k], implicitarcs)

            if !isempty(implicitarcs)
                # If there are implicit arcs in the optimal dual path, add the corresponding dual
                dualoptimals[k] = false
                union!(explicitarcs[k], implicitarcs)
                MOI.submit.(model, MOI.LazyConstraint(cb_data), dualfeas[k][implicitarcs])
            end
        end

        # Print
        for k in 1:length(prob.K)
            primalcolor = primaloptimals[k] ? GREEN_FG : YELLOW_FG
            print(primalcolor, @sprintf("%2d", length(explicitnodes[k])))
            dualcolor = dualoptimals[k] ? GREEN_FG : YELLOW_FG
            print(dualcolor, @sprintf("(%2d) ", length(explicitarcs[k])))
        end
        println(crayon"reset")
    end
end
