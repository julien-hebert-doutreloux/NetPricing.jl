## Build graph from problem
function build_graph(prob::AbstractProblem)
    g = SimpleWeightedDiGraph(nodes(prob))
    reset!(g, prob)
    return g
end

## Enable/disable arcs/nodes
function reset!(graph::AbstractGraph, prob::AbstractProblem)
    for arc in arcs(prob)
        add_edge!(graph, arc.src, arc.dst, arc.cost + eps(0.0))
    end
    return graph
end

function restore!(graph::AbstractGraph, weights)
    @inbounds graph.weights = copy(weights)
end

function disable_arcs!(graph::AbstractGraph, prob::AbstractProblem, _arcs)
    for arc in @view arcs(prob)[_arcs]
        @inbounds graph.weights[arc.dst, arc.src] = Inf
    end
end

function disable_nodes!(graph::AbstractGraph, nodes)
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
    if orig == dest
        return Int[], 0.0
    end

    dists = fill(Inf, nv(graph))
    parents = zeros(Int, nv(graph))

    ws = Graphs.weights(graph)
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

# Path arc iterator
struct PathArcIterator
    path::Vector{Int}
    arcdict::Dict{Tuple{Int,Int},Int}
end

path_arcs(path, arcdict::Dict{Tuple{Int,Int},Int}) = PathArcIterator(path, arcdict)
path_arcs(path, prob::AbstractProblem) = PathArcIterator(path, srcdst_to_index(prob))

function Base.iterate(iter::PathArcIterator, i = 0)
    @unpack path, arcdict = iter
    return length(path) - i > 1 ? (arcdict[(path[i + 1],path[i + 2])], i + 1) : nothing
end

Base.eltype(::Type{PathArcIterator}) = Int
Base.length(iter::PathArcIterator) = length(iter.path) - 1

# Get path cost
function get_path_cost(path, arccosts)
    cost = 0.0
    for i in 1:(length(path)-1)
        cost += arccosts[(path[i], path[i+1])]
    end
    return cost
end
get_path_cost(path, prob::AbstractProblem) = get_path_cost(path, srcdst_to_cost(prob))
