function preprocess_light(prob::Problem, k)
    graph = build_graph(prob)
    revgraph = reverse(graph)
    infgraph = copy(graph)
    disable_tolls!(infgraph, prob)

    nv = prob.V

    orig, dest = prob.K[k].orig, prob.K[k].dest
    omin = dijkstra_shortest_paths(graph, orig).dists
    dmin = dijkstra_shortest_paths(revgraph, dest).dists
    _, odmax = shortest_path(infgraph, orig, dest)

    relevants = BitSet(i for i in 1:nv if omin[i] + dmin[i] <= odmax || i == orig || i == dest)

    Vmap = collect(relevants)
    Amap = findall(arc -> arc.src ∈ relevants && arc.dst ∈ relevants, prob.A)

    Vrevmap = revmap(Vmap, prob.V)
    Arevmap = revmap(Amap, length(prob.A))

    # Translate the arcs
    A = [ProblemArc(Vrevmap[arc.src], Vrevmap[arc.dst], arc.cost, arc.toll) for arc in @view prob.A[Amap]]

    return PreprocessedProblem(prob, length(Vmap), A, Vmap, Amap, Vrevmap, Arevmap, k,
        Vrevmap[prob.K[k].orig],
        Vrevmap[prob.K[k].dest])
end
