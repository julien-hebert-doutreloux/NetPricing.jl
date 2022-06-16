function generate_problem(graph::SimpleWeightedDiGraph, num_commodities, args::ProblemGenerationArgs)
    (; demand_dist) = args
    
    # Randomize commodities
    V = nv(graph)
    odpairs = sample([(j, i) for i in 1:V, j in 1:V if i != j], num_commodities, replace=false)
    K = [Commodity(odpairs[k]..., rand(demand_dist)) for k in 1:num_commodities]

    # Generate tolled arcs
    A = generate_tolled_arcs(graph, K, args)

    return Problem(V, A, K)
end

generate_problem(graph::AbstractGraph, num_commodities, args) =
    generate_problem(generate_cost(graph, args), num_commodities, args)

generate_problem(graph, num_commodities; kwargs...) = generate_problem(graph, num_commodities, ProblemGenerationArgs(; kwargs...))

function generate_cost(graph::SimpleDiGraph, args::ProblemGenerationArgs)
    (; cost_dist, max_cost, max_cost_prob, symmetric_cost) = args
    
    weighted_graph = SimpleWeightedDiGraph(nv(graph))
    for edge in edges(graph)
        i, j = src(edge), dst(edge)
        cost = _random_arc_cost(max_cost_prob, max_cost, cost_dist)

        add_edge!(weighted_graph, i, j, cost + eps(0.0))
        symmetric_cost && has_edge(weighted_graph, j, i) && add_edge!(weighted_graph, j, i, cost + eps(0.0))
    end
    return weighted_graph
end
generate_cost(graph::AbstractGraph, args) = generate_cost(SimpleDiGraph(graph), args)
generate_cost(graph; kwargs...) = generate_cost(graph, ProblemGenerationArgs(; kwargs...))

_random_arc_cost(max_cost_prob, max_cost, cost_dist) = rand() < max_cost_prob ? max_cost : rand(cost_dist)

function generate_tolled_arcs(graph::SimpleWeightedDiGraph, commodities, args::ProblemGenerationArgs)
    (; tolled_proportion, symmetric_tolled, random_tolled_proportion) = args
    
    # Shortest path statistics
    occurences = Dict((src(e), dst(e)) => 0 for e in edges(graph))
    for comm in commodities
        path, _  = shortest_path(graph, comm.orig, comm.dest) 
        for pair in consecutive_pairs(path)
            occurences[pair] += 1
        end
    end

    # Merge with reversed arc if symmetric_tolled = true
    if symmetric_tolled
        for ((i, j), _) in occurences
            (i > j) && continue                         # Only do once per pair
            haskey(occurences, (j, i)) || continue      # Must have reversed arc
            occurences[(i, j)] = occurences[(j, i)] = occurences[(i, j)] + occurences[(j, i)]
        end
    end

    # Sort the occurences (reversed order)
    candidates = first.(collect(sort(occurences, by=k->occurences[k])))

    # Toll-free graph to check if all commodities has at least 1 toll-free path
    tollfreegraph = SimpleDiGraph(graph)

    tolled = Set{Tuple{Int,Int}}()
    num_arcs = ne(graph)
    num_tolled_arcs = round(Int, num_arcs * tolled_proportion)
    num_tolled_from_sorted = round(Int, num_tolled_arcs * (1 - random_tolled_proportion))
    shuffled = false

    # Convert arcs to tolled arcs
    while !isempty(candidates)
        i, j = pop!(candidates)
        (i, j) in tolled && continue                            # Already in tolled

        _is_removable(tollfreegraph, i, j, commodities) || continue        # Must be removable (does not disconnect any commodity)
        if symmetric_tolled && has_edge(tollfreegraph, j, i)
            _is_removable(tollfreegraph, j, i, commodities) || continue    # So does the reversed arc if symmetric_tolled
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

        # Sorted target reached: shuffle the remaining candidates
        if !shuffled && length(tolled) >= num_tolled_from_sorted
            shuffle!(candidates)
            shuffled = true
        end
    end

    # Convert to ProblemArc
    A = _convert_to_A(graph, tolled)

    return A
end
generate_tolled_arcs(graph, commodities; kwargs...) = generate_tolled_arcs(graph, commodities, ProblemGenerationArgs(; kwargs...))

function _convert_to_A(graph, tolled)
    A = ProblemArc[]
    for edge in edges(graph)
        i, j = src(edge), dst(edge)
        cost = edge.weight
        toll = (i, j) in tolled
        toll && (cost /= 2)
        push!(A, ProblemArc(i, j, cost, toll))
    end
    return A
end

function _is_removable(graph, i, j, commodities)
    rem_edge!(graph, i, j)
    removable = all(comm->has_path(graph, comm.orig, comm.dest), commodities)
    add_edge!(graph, i, j)
    return removable
end