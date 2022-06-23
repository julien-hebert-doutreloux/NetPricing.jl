# Generate a list of problems which are highly correlated to each other from max_size to min_size.
# The next problem reuses the grid, the set of tolled arcs and a modified set of commodities of the previous problem.
# This procedure is used to generate problems to evaluate the effects of strong bilevel feasibility cuts.
function generate_progressive_grid(max_size, min_size, base_num_commodities, args::ProblemGenerationArgs)
    # Generate the base problem (max_size)
    graph = generate_cost(grid_graph(max_size, max_size), args)

    num_commodities = base_num_commodities
    prob = generate_problem(graph, num_commodities, args)
    problems = [prob]

    # Generate the smaller problems progressively
    for size in (max_size-1):-1:min_size
        # Randomly choose the corner to reuse
        offset_x, offset_y = rand(0:1, 2)   # offset_x means left/right, offset_y means top/bottom

        # Remap the nodes
        node_map, inv_node_map = Dict(), Dict()
        grid_index(x, y, s) = (y - 1) * s + x
        for x in 1:size, y in 1:size
            prev_ind = grid_index(x + offset_x, y + offset_y, size + 1)
            next_ind = grid_index(x, y, size)
            node_map[prev_ind] = next_ind
            inv_node_map[next_ind] = prev_ind
        end
        subgraph_nodes = BitSet(keys(node_map))

        # Extract the induced subgraph (with the old costs of arcs)
        next_graph = graph[collect(subgraph_nodes)]
        V = nv(next_graph)

        # Increase the number of commodities inversely-proportionally to the number of edges
        # The number of edges is 2n(n-1) where n is the size of the grid
        num_commodities = round(Int, base_num_commodities * (max_size * (max_size - 1)) / (size * (size - 1)))
        K = _generate_progressive_commodities(graph, next_graph, prob.K, node_map, num_commodities, args)

        # Add the tolled arcs (from scratch)
        A = generate_tolled_arcs(next_graph, K, args)

        # Update variables for the next iterations
        prob = Problem(V, A, K)
        push!(problems, prob)
        graph = next_graph
    end

    return problems
end

generate_progressive_grid(max_size, min_size, base_num_commodities; kwargs...) =
    generate_progressive_grid(max_size, min_size, base_num_commodities, ProblemGenerationArgs(; kwargs...))

function generate_progressive_delaunay(num_nodes_array, base_num_commodities, args::ProblemGenerationArgs)
    # Generate the base problem
    raw_graph, seeds = _delaunay_graph_with_coors(num_nodes_array[1])
    graph = generate_cost(raw_graph, args)

    num_commodities = base_num_commodities
    prob = generate_problem(graph, num_commodities, args)
    problems = [prob]

    # Routine to calculate the central seed (used as the first node in the augmentation)
    function get_central_seed()
        avg_x = sum(getx, seeds) / length(seeds)
        avg_y = sum(gety, seeds) / length(seeds)
        return argmin(map(p -> (getx(p) - avg_x)^2 + (gety(p) - avg_y)^2, seeds))
    end
    central = get_central_seed()

    # Generate the smaller problems progressively
    for num_nodes in num_nodes_array[2:end]
        # Remap the nodes
        next_nodes = _augment_delaunay(graph, num_nodes, central)
        node_map = Dict(reverse.(enumerate(next_nodes)))

        # Make a new graph from the subset of seeds, reuse the old set of arc costs
        next_graph = graph[next_nodes]
        V = nv(next_graph)

        # Increase the number of commodities inversely-proportionally to the number of edges
        # For Delaunay graph, we don't know the number of edges from the number of nodes
        num_commodities = round(Int, base_num_commodities * ne(raw_graph) * 2 / ne(next_graph))
        K = _generate_progressive_commodities(graph, next_graph, prob.K, node_map, num_commodities, args)

        # Add the tolled arcs (from scratch)
        A = generate_tolled_arcs(next_graph, K, args)

        # Update variables for the next iterations
        prob = Problem(V, A, K)
        push!(problems, prob)
        graph = next_graph

        # Recalculate central seed
        seeds = seeds[next_nodes]
        central = get_central_seed()
    end

    return problems
end

generate_progressive_delaunay(num_nodes_array, base_num_commodities; kwargs...) =
    generate_progressive_delaunay(num_nodes_array, base_num_commodities, ProblemGenerationArgs(; kwargs...))

# Sub-routines
# Reuse old set of commodities
function _generate_progressive_commodities(last_graph, next_graph, last_K, node_map, num_commodities, args::ProblemGenerationArgs)
    (; demand_dist) = args

    subgraph_nodes = BitSet(keys(node_map))
    next_K = Commodity[]

    # Reuse the commodities, if the origin/destination is outside of the new graph, find the outer-most node on the shortest path.
    V = nv(next_graph)
    available_odpairs = Set([(o, d) for o in 1:V, d in 1:V if o != d])  # Keep track of valid od-pairs
    for comm in last_K
        o, d  = comm.orig, comm.dest
        # If origin or destination are not in the subgraph, find the shortest path,
        # then use the first and last nodes in the subgraph as origin and destination
        if o ∉ subgraph_nodes || d ∉ subgraph_nodes
            path, _ = shortest_path(last_graph, o, d)
            oi = findfirst(in(subgraph_nodes), path)
            di = findlast(in(subgraph_nodes), path)
            # If the shortest path is not in the subgraph, skip this commodity
            if isnothing(oi) || isnothing(di)
                continue
            end
            o, d = path[oi], path[di]
        end
        # Remap (o, d) to the new graph
        o, d = node_map[o], node_map[d]
        # Check if (o, d) is valid
        ((o, d) in available_odpairs) || continue
        # Otherwise, add to the new commodity list
        push!(next_K, Commodity(o, d, comm.demand))
        delete!(available_odpairs, (o, d))
    end

    # Add new commodities to reach the required number of commodities
    sampled_odpairs = sample(collect(available_odpairs), num_commodities - length(next_K), replace=false)
    for (o, d) in sampled_odpairs
        push!(next_K, Commodity(o, d, rand(demand_dist)))
    end

    return next_K
end

# Reuse old set of tolled arcs (obsolete)
function _generate_progressive_tolled_arcs(graph, last_A, node_map, args::ProblemGenerationArgs)
    (; tolled_proportion, symmetric_tolled) = args
    subgraph_nodes = BitSet(keys(node_map))

    # Try to reuse the tolled arcs
    num_arcs = ne(graph)
    num_tolled_arcs = round(Int, num_arcs * tolled_proportion)

    candidates = [(arc.src, arc.dst) for arc in last_A if arc.toll] # Extract tolled arcs from previous problem
    filter!(candidates) do (i, j)                                   # Filter out the arcs that are not in the subgraph
        return i in subgraph_nodes && j in subgraph_nodes
    end
    map!(candidates, candidates) do (i, j)                          # Map to subgraph's node indices
        return (node_map[i], node_map[j])
    end
    shuffle!(candidates)

    # Append the remaining arcs to candidates in random order
    remaining = [(src(e), dst(e)) for e in edges(graph)]
    setdiff!(remaining, candidates)
    shuffle!(remaining)
    append!(candidates, remaining)

    # Add the candidates (old tolled arcs first, then other random arcs) if they do not disconnect any o-d pair
    tollfreegraph = SimpleDiGraph(graph)
    tolled = Set{Tuple{Int,Int}}()

    reverse!(candidates)    # Reverse the order so we can use pop!
    while !isempty(candidates)
        i, j = pop!(candidates)
        # Already in tolled
        (i, j) in tolled && continue

        # Must be removable (does not disconnect any commodity)
        _is_removable(tollfreegraph, i, j, next_K) || continue
        # So does the reversed arc if symmetric_tolled
        if symmetric_tolled && has_edge(tollfreegraph, j, i)
            _is_removable(tollfreegraph, j, i, next_K) || continue
        end

        # If OK, add to tolled
        push!(tolled, (i, j))
        rem_edge!(tollfreegraph, i, j)
        if symmetric_tolled && has_edge(tollfreegraph, j, i)
            push!(tolled, (j, i))
            rem_edge!(tollfreegraph, j, i)
        end

        # Total target reached: break
        length(tolled) >= num_tolled_arcs && break
    end

    A = _convert_to_A(graph, tolled)
    return A
end

# Remove some nodes from a Delaunay graph, returns the list of the remaining nodes (obsolete)
function _reduce_delaunay(graph, target)
    graph = copy(graph)
    remaining_nodes = collect(1:nv(graph))
    while nv(graph) > target
        # Remove the node with the smallest degree
        i = argmin(degree(graph))
        rem_vertex!(graph, i)
        deleteat!(remaining_nodes, i)
    end
    return remaining_nodes
end

# Fill in the old costs or generate new one for new pairs (obsolete)
function _fill_costs(graph::SimpleDiGraph, arc_costs_dict, args::ProblemGenerationArgs)
    (; symmetric_cost) = args
    
    weighted_graph = SimpleWeightedDiGraph(nv(graph))
    for edge in edges(graph)
        i, j = src(edge), dst(edge)
        if !haskey(arc_costs_dict, (i, j))
            c = _random_arc_cost(args)
            arc_costs_dict[(i, j)] = c
            symmetric_cost && (arc_costs_dict[(j, i)] = c)
        end
        cost = arc_costs_dict[(i, j)]
        add_edge!(weighted_graph, i, j, cost + eps(0.0))
        symmetric_cost && has_edge(weighted_graph, j, i) && add_edge!(weighted_graph, j, i, cost + eps(0.0))
    end
    return weighted_graph
end
_fill_costs(graph::AbstractGraph, arc_costs_dict, args) = _fill_costs(SimpleDiGraph(graph), arc_costs_dict, args)

# Augment Delaunay graph
function _augment_delaunay(graph, target, first)
    chosen = [first]
    remaining = setdiff(1:nv(graph), first)
    for _ in 2:target
        # We add the node with the highest number of connections to the current subset of nodes
        # If there're ties, choose the node with highest degree
        criteria = map(remaining) do i
            num_connections = count(j -> has_edge(graph, i, j), chosen)     # Number of connections to to current subset
            node_degree = degree(graph, i)
            return num_connections, node_degree
        end

        next_node = remaining[argmax(criteria)]
        push!(chosen, next_node)
        setdiff!(remaining, next_node)
    end
    return chosen
end
