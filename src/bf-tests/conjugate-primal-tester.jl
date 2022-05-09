struct ConjugatePrimalTester <: AbstractBFTester
    model::Model
    prob::Problem
    num_commodities::Int
end

# Interface
problem(tester::ConjugatePrimalTester) = tester.prob
num_commodities(tester::ConjugatePrimalTester) = tester.num_commodities

# Init a conjugate model to test bilevel feasibility between some paths, weights = 1
function ConjugatePrimalTester(prob::Problem, num_commodities;
    silent=true, threads=nothing, sdtol=default_sdtol())

    model = Model(() -> Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), silent)
    set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)

    nk = num_commodities
    na = length(arcs(prob))
    a1 = tolled_arcs(prob)

    c = cost_vector(prob)
    A = incidence_matrix(prob)

    @variable(model, 0 ≤ x[i=1:na, k=1:nk] ≤ 1)
    @variable(model, L)
    @objective(model, Max, 0)

    @constraint(model, primalfeas, A * x .== 0)
    @constraint(model, demandlimit[a=a1], sum(x[a,:]) ≤ 0)
    @constraint(model, strongdual, sum(c' * x) ≤ L + sdtol)
    
    return ConjugatePrimalTester(model, prob, num_commodities)
end

# Set OD pairs
function set_odpairs(tester::ConjugatePrimalTester, odpairs::AbstractVector{Tuple{Int,Int}})
    model, prob = tester.model, tester.prob
    b = hcat([sourcesink_vector(prob, o, d) for (o, d) in odpairs]...)
    set_normalized_rhs.(model[:primalfeas], b)
    return
end

# Set paths
function set_paths(tester::ConjugatePrimalTester, paths; set_odpairs=true)
    prob, model = tester.prob, tester.model
    set_odpairs && NetPricing.set_odpairs(tester, [(path[1], path[end]) for path in paths])
    counts = count_tolled_arcs(prob, paths)
    for a in tolled_arcs(prob)
        set_normalized_rhs(model[:demandlimit][a], get(counts, a, 0.))
    end
end

# Bilevel feasibility test
function is_bilevel_feasible(tester::ConjugatePrimalTester, paths; set_odpairs=true)
    prob, model = tester.prob, tester.model
    set_paths(tester, paths, set_odpairs=set_odpairs)

    # Find the minimum cost using L
    L = model[:L]
    is_fixed(L) && unfix(L)
    @objective(model, Min, L)

    # Solve
    optimize!(model)

    # Compare with input costs
    objval = objective_value(model)
    expecting = sum_paths_costs(prob, paths)
    
    return objval >= expecting - 1e-6
end

# Strong bilevel feasibility test
function is_strongly_bilevel_feasible(tester::ConjugatePrimalTester, paths; set_odpairs=true)
    prob, model = tester.prob, tester.model
    nk = tester.num_commodities
    na = length(arcs(prob))

    set_paths(tester, paths, set_odpairs=set_odpairs)
    
    Lval = sum_paths_costs(prob, paths)
    fix(model[:L], Lval, force=true)
    
    # Extract the set of arcs for each path
    arcdict = srcdst_to_index(prob)
    active_arcs = [BitSet(path_arcs(path, arcdict)) for path in paths]

    # Maximize x of arcs that do not appear in active_arcs
    x = model[:x]
    @objective(model, Max, sum(x[a,k] for k in 1:nk, a in 1:na if a ∉ active_arcs[k]))

    # Solve
    # If obj = 0, the paths are strongly bilevel feasible
    # If obj > 0, they are either bilevel infeasible or weakly bilevel feasible
    optimize!(model)
    return objective_value(model) ≤ 1e-6
end
