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
function set_tolls!(graph::AbstractGraph, preprob::PreprocessedProblem, tolls::AbstractTollsDict)
    # Translate tolls dict to preprob
    Arevmap = preprob.Arevmap
    mappedtolls = Dict(Arevmap[index] => toll for (index, toll) in tolls if Arevmap[index] > 0)
    return _set_tolls_abstract!(graph, preprob, mappedtolls)
end

set_tolls!(graph::AbstractGraph, preprob::PathPreprocessedProblem, tolls::AbstractTollsDict) = set_tolls!(graph, preprob.pprob, tolls)

## Reset tolls to 0
reset_tolls!(graph::AbstractGraph, prob::AbstractProblem) = _set_tolls_abstract!(graph, prob, zero_tolls(prob))

## Disable all tolled arcs
disable_tolls!(graph::AbstractGraph, prob::AbstractProblem) = _set_tolls_abstract!(graph, prob, inf_tolls(prob))