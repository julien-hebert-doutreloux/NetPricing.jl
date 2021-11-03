## Build graph from problem
function build_graph(prob::Problem)
    g = SimpleWeightedDiGraph(prob.V)
    reset!(g, prob)
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

# Enable/disable arcs/nodes
function reset!(graph::AbstractGraph, prob::Problem)
    for arc in prob.A
        add_edge!(graph, arc.src, arc.dst, arc.cost + eps(0.0))
    end
    return graph
end

function restore!(graph::AbstractGraph, weights)
    @inbounds graph.weights = copy(weights)
end

function disable_arcs!(graph::AbstractGraph, prob::Problem, arcs::Vector{Int})
    for arc in @view prob.A[arcs]
        @inbounds graph.weights[arc.dst, arc.src] = Inf
    end
end

function disable_nodes!(graph::AbstractGraph, nodes::Vector{Int})
    for i in nodes
        for j in outneighbors(graph, i)
            @inbounds graph.weights[j, i] = Inf
        end
    end
end

## Find the shortest path
function shortest_path_old(graph::AbstractGraph, orig, dest)
    ds = dijkstra_shortest_paths(graph, orig)
    path = enumerate_paths(ds, dest)
    return path, ds.dists[dest]
end

function shortest_path(graph::AbstractGraph, orig, dest)
    dists = fill(Inf, nv(graph))
    parents = zeros(Int, nv(graph))

    ws = weights(graph)
    queue = PriorityQueue(orig => 0.0)
    dists[orig] = 0

    # Forward
    while !isempty(queue)
        u = dequeue!(queue)
        (u == dest) && break

        dist_u = dists[u]

        for v in outneighbors(graph, u)
            newdist = dist_u + ws[u,v]
            if newdist < dists[v]
                queue[v] = dists[v] = newdist
                parents[v] = u
            end
        end
    end

    # Backward
    (parents[dest] == 0) && return Int[], Inf

    path = [dest]
    v = dest
    while v != orig
        v = parents[v]
        push!(path, v)
    end

    return reverse!(path), dists[dest]
end
