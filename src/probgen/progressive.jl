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
    base_graph = generate_cost(raw_graph, args)

    # Function to evaluate number of nodes inside a window
    max_width = max_coord - min_coord
    mid_coord = (max_coord + min_coord) / 2
    in_window(x::Number, width) = mid_coord - width / 2 <= x <= mid_coord + width / 2
    in_window(p::Point2D, width) = in_window(getx(p), width) && in_window(gety(p), width)
    count_nodes_window(width) = count(p -> in_window(p, width), seeds)

    num_commodities = base_num_commodities
    prob = generate_problem(base_graph, num_commodities, args)
    problems = [prob]
    last_inv_node_map = Dict(i => i for i in 1:num_nodes_array[1])

    # Generate the smaller problems progressively
    for num_nodes in num_nodes_array[2:end]
        # Find the window that contains the target number of nodes
        count_fn = w -> count_nodes_window(w) - num_nodes
        window = binary_root_search(count_fn, 0, max_width)

        # Remap the nodes
        next_nodes = [i for (i, p) in enumerate(seeds) if in_window(p, window)]
        node_map = Dict(reverse.(enumerate(next_nodes)))

        # Extract the induced subgraph (with the old costs of arcs)
        next_graph = base_graph[next_nodes]
        V = nv(next_graph)

        # Increase the number of commodities inversely-proportionally to the number of edges
        # For Delaunay graph, we don't know the number of edges from the number of nodes
        num_commodities = round(Int, base_num_commodities * ne(base_graph) / ne(next_graph))
        # We need to translate the last set of commodities to the indices of the base graph
        last_K = [Commodity(last_inv_node_map[comm.orig], last_inv_node_map[comm.dest], comm.demand) for comm in prob.K]
        K = _generate_progressive_commodities(base_graph, next_graph, last_K, node_map, num_commodities, args)

        # Add the tolled arcs (from scratch)
        A = generate_tolled_arcs(next_graph, K, args)

        # Update variables for the next iterations
        prob = Problem(V, A, K)
        push!(problems, prob)
        last_inv_node_map = Dict(enumerate(next_nodes))
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