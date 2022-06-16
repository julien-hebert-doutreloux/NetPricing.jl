# Generate a list of problems which are highly correlated to each other from max_size to min_size.
# The next problem reuses the grid, the set of tolled arcs and a modified set of commodities of the previous problem.
# This procedure is used to generate problems to evaluate the effects of strong bilevel feasibility cuts.
function generate_progressive_grid(max_size, min_size, base_num_commodities, args::ProblemGenerationArgs)
    (; demand_dist, tolled_proportion, symmetric_tolled) = args

    # Generate the base problem (max_size)
    graph = generate_cost(grid_graph(max_size, max_size), args)

    num_commodities = base_num_commodities
    prob = generate_problem(graph, num_commodities, args)
    problems = [prob]

    # Generate the smaller problems prograssively
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

        # Increase the number of commodities inversely-proportionally to the number of nodes
        num_commodities = round(Int, base_num_commodities * (max_size * max_size) / (size * size))
        last_K = prob.K
        next_K = Commodity[]

        # Reuse the commodities, if the origin/destination is outside of the new graph, find the outer-most node on the shortest path.
        V = nv(next_graph)
        available_odpairs = Set([(o, d) for o in 1:V, d in 1:V if o != d])  # Keep track of valid od-pairs
        for comm in last_K
            o, d  = comm.orig, comm.dest
            # If origin or destination are not in the subgraph, find the shortest path,
            # then use the first and last nodes in the subgraph as origin and destination
            if o ∉ subgraph_nodes || d ∉ subgraph_nodes
                path, _ = shortest_path(graph, o, d)
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

        # Try to reuse the tolled arcs
        num_arcs = ne(next_graph)
        num_tolled_arcs = round(Int, num_arcs * tolled_proportion)

        candidates = [(arc.src, arc.dst) for arc in prob.A if arc.toll] # Extract tolled arcs from previous problem
        filter!(candidates) do (i, j)                                   # Filter out the arcs that are not in the subgraph
            return i in subgraph_nodes && j in subgraph_nodes
        end
        map!(candidates, candidates) do (i, j)                          # Map to subgraph's node indices
            return (node_map[i], node_map[j])
        end
        shuffle!(candidates)

        # Append the remaining arcs to candidates in random order
        remaining = [(src(e), dst(e)) for e in edges(next_graph)]
        setdiff!(remaining, candidates)
        shuffle!(remaining)
        append!(candidates, remaining)

        # Add the candidates (old tolled arcs first, then other random arcs) if they do not disconnect any o-d pair
        tollfreegraph = SimpleDiGraph(next_graph)
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

        # Add the tolled arcs (from scratch)
        A = _convert_to_A(next_graph, tolled)

        # Update variables for the next iterations
        prob = Problem(V, A, next_K)
        push!(problems, prob)
        graph = next_graph
    end

    return problems
end

generate_progressive_grid(max_size, min_size, base_num_commodities; kwargs...) =
    generate_progressive_grid(max_size, min_size, base_num_commodities, ProblemGenerationArgs(; kwargs...))