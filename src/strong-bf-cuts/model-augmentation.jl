# Add strong bf cuts for some pairs of formulations
function add_strong_bf_cuts(model::Model, forms::Vector{<:Formulation}; kwargs...)
    primals = convert.(PrimalPath, filter!(p -> p isa PrimalPath, primal.(forms)))
    add_strong_bf_cuts(model, primals; kwargs...)
end

function add_strong_bf_cuts(model::Model, primals::Vector{PrimalPath}; maxpairs=1000, commpairs=100, kwargs...)
    all_commpairs = collect(subsets(primals, 2))
    filter!(p -> length(p[1].prob.paths) * length(p[2].prob.paths) <= maxpairs, all_commpairs)
    sort!(all_commpairs, by=p->commpair_score(p[1].prob, p[2].prob), rev=true)

    cuts = ConstraintRef[]
    for (p1, p2) in all_commpairs[1:min(end, commpairs)]
        append!(cuts, generate_strong_bf_cuts(model, p1, p2; kwargs...))
    end
    return cuts
end

# Score to rank the pairs of commodities
function commpair_score(prob1::PathPreprocessedProblem, prob2::PathPreprocessedProblem)
    # Count the number of tolled arcs in common
    prob1_tolled = BitSet(arcmap(prob1)[tolled_arcs(prob1)])
    prob2_tolled = BitSet(arcmap(prob2)[tolled_arcs(prob2)])
    numcommon = length(prob1_tolled ∩ prob2_tolled)

    # Discourage prob with many solutions
    return numcommon / log(length(prob1.paths) * length(prob2.paths))^2
end
