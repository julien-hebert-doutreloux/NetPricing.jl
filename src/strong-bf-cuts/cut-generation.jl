# Generate strong bilevel feasibility cuts: test all pairs of paths of 2 commodities,
# then find a set of cuts that covers the non-strongly-bf combinations
function generate_strong_bf_cuts(model::Model, primal1::PrimalPath, primal2::PrimalPath; kwargs...)
    sbfmap = strong_bf_map(primal1.prob, primal2.prob; kwargs...)     # Make the strong-bf map of all pairs
    covering = biclique_cover(.!sbfmap)                  # Find a cover of all non-strong-bf pairs

    # Convert the covering to cuts
    z1, z2 = primal1.z, primal2.z
    all_constraints = ConstraintRef[]
    for (l, r) in covering
        push!(all_constraints, @constraint(model, sum(z1[l]) + sum(z2[r]) <= 1))
    end
    return all_constraints
end

# Make a BitMatrix indicating which pairs of paths are strongly-bf
function strong_bf_map(prob1::PathPreprocessedProblem, prob2::PathPreprocessedProblem; threads=nothing, bf_tester=ConjugatePrimalTester)
    np, nq = length(prob1.paths), length(prob2.paths)
    sbfmap = BitMatrix(undef, np, nq)

    tester = bf_tester(parent(prob1), 2; threads=threads)
    set_odpairs(tester, [index(prob1), index(prob2)])

    for (i, p) in enumerate(paths(prob1)), (j, q) in enumerate(paths(prob2))
        sbfmap[i, j] = is_strongly_bilevel_feasible(tester, [p, q], set_odpairs=false)
    end

    return sbfmap
end