function _set_tolls_abstract!(graph::AbstractGraph, prob::AbstractProblem, tolls::AbstractTollsDict)
    for (index, toll) in tolls
        arc = arcs(prob)[index]
        arc.toll || throw(ErrorException("Arc $index is not a tolled arc"))

        add_edge!(graph, arc.src, arc.dst, arc.cost + toll + eps(0.0))
    end
    return graph
end

## Set tolls to unprocessed graph
set_tolls!(graph::AbstractGraph, prob::Problem, tolls::AbstractTollsDict) = _set_tolls_abstract!(graph, prob, tolls)

## Set tolls to preprocessed graph
function set_tolls!(graph::AbstractGraph, prob::PreprocessedProblem, tolls::AbstractTollsDict)
    # Translate tolls dict to prob
    Arevmap = prob.Arevmap
    for (index, toll) in tolls
        mapped_index = Arevmap[index]
        iszero(mapped_index) && continue

        arc = arcs(prob)[mapped_index]
        arc.toll || throw(ErrorException("Arc $mapped_index is not a tolled arc"))
        add_edge!(graph, arc.src, arc.dst, arc.cost + toll + eps(0.0))
    end
    return graph
end

set_tolls!(graph::AbstractGraph, prob::PathPreprocessedProblem, tolls::AbstractTollsDict) = set_tolls!(graph, prob.pprob, tolls)

## Reset tolls to 0
reset_tolls!(graph::AbstractGraph, prob::AbstractProblem) = _set_tolls_abstract!(graph, prob, zero_tolls(prob))

## Disable all tolled arcs
disable_tolls!(graph::AbstractGraph, prob::AbstractProblem) = _set_tolls_abstract!(graph, prob, inf_tolls(prob))