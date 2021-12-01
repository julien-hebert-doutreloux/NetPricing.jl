function adaptive_info(prob::Problem)
    graph = build_graph(prob)
    infgraph = copy(graph)
    disable_tolls!(infgraph, prob)

    nk = length(prob.K)
    nv = prob.V

    all_relevants = BitSet[]
    omin = zeros(nk, nv)
    dmin = zeros(nk, nv)
    all_odmax = zeros(nk)   # Min toll-free dist from o to d

    for (k, comm) in enumerate(prob.K)
        orig, dest = comm.orig, comm.dest
        _, odmax = shortest_path(infgraph, orig, dest)

        revgraph = reverse(graph)
        omin[k, :] = dijkstra_shortest_paths(graph, orig).dists
        dmin[k, :] = dijkstra_shortest_paths(revgraph, dest).dists
    
        relevants = BitSet(i for i in 1:nv if omin[k,i] + dmin[k,i] <= odmax || i == orig || i == dest)
        push!(all_relevants, relevants)
        all_odmax[k] = odmax
    end

    return all_relevants, omin, dmin, all_odmax
end

function calculate_odmax(prob::Problem)
    graph = build_graph(prob)
    disable_tolls!(graph, prob)

    odmax = zeros(length(prob.K))
    for (k, comm) in enumerate(prob.K)
        orig, dest = comm.orig, comm.dest
        _, odmax[k] = shortest_path(graph, orig, dest)
    end

    return odmax
end

function calculate_mincost_fixedarc(prob::Problem)
    graph = build_graph(prob)
    revgraph = reverse(graph)

    mincost = zeros(length(prob.K), length(prob.A))
    for (k, comm) in enumerate(prob.K)
        orig, dest = comm.orig, comm.dest

        omin = dijkstra_shortest_paths(graph, orig).dists
        dmin = dijkstra_shortest_paths(revgraph, dest).dists

        for a in 1:length(prob.A)
            arc = prob.A[a]
            mincost[k, a] = omin[arc.src] + arc.cost + dmin[arc.dst]
        end
    end

    return mincost
end
