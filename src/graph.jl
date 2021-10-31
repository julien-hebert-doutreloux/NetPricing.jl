## Build graph from problem
function build_graph(prob::Problem)
    g = SimpleWeightedDiGraph(prob.V)
    for arc in prob.A
        add_edge!(g, arc.src, arc.dst, arc.cost + eps(0.0))
    end
    return g
end

## Set tolls to graph
function set_tolls!(graph::AbstractGraph, prob::Problem, tolls::AbstractTollsDict)
    for (index, toll) in tolls
        arc = prob.A[index]
        arc.toll || throw(ErrorException("Arc $index is not a tolled arc"))

        add_edge!(graph, arc.src, arc.dst, arc.cost + toll + eps(0.0))
    end
    return graph
end

reset_tolls!(graph::AbstractGraph, prob::Problem) = set_tolls!(graph, prob, zero_tolls(prob))

disable_arcs!(graph::AbstractGraph, prob::Problem, arcs) = set_tolls!(graph, prob, Dict(arcs .=> Inf))

## Find the shortest path
function shortest_path(graph::AbstractGraph, orig, dest)
    ds = dijkstra_shortest_paths(graph, orig)
    path = enumerate_paths(ds, dest)
    return path, ds.dists[dest]
end
