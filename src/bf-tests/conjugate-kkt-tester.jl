# ConjugateKKTModel does support strong bilevel feasibility test
const ConjugateKKTTester = ConjugateBFTester{ConjugateKKTModel}

# Strong bilevel feasible test
function is_strongly_bilevel_feasible(tester::ConjugateKKTTester, paths; set_odpairs=true)
    cmodel = tester.solver
    prob, model = problem(cmodel), cmodel.model

    # Set demands from paths, using SlackPriority
    set_odpairs && NetPricing.set_odpairs(tester, [(path[1], path[end]) for path in paths])
    set_demands(cmodel, count_tolled_arcs(prob, paths), SlackPriority())
    
    # Extract the set of arcs for each path
    arcdict = srcdst_to_index(prob)
    active_arcs = [BitSet(path_arcs(path, arcdict)) for path in paths]

    # Set the model into interior test mode
    set_slack(cmodel, active_arcs)

    # Solve
    optimize!(model)

    # Strong bilevel feasibility implies S > 0
    # If S < 0, the paths are not bilevel feasible
    # If S = 0, they are bilevel feasible but not strongly
    Sval = objective_value(model)
    return Sval ≥ 1e-6
end
