function calculate_bigM(prob::AbstractCommodityProblem)
    parentprob = parent(prob)
    graph = build_graph(parentprob)

    a1 = tolled_arcs(parentprob)
    comm = parentprob.K[index(prob)]
    Amap = BitSet(arcmap(prob))

    reset_tolls!(graph, parentprob)
    omin = dists_from(graph, comm.orig)
    dmin = dists_to(graph, comm.dest)

    disable_tolls!(graph, parentprob)
    omax = dists_from(graph, comm.orig)
    dmax = dists_to(graph, comm.dest)
    odmax = omax[comm.dest]

    p = [shortest_path(graph, a.src, a.dst)[2] for a in parentprob.A[a1]]

    M = zeros(length(a1))
    for (ia, a) in enumerate(a1)
        (a ∉ Amap) && continue      # Unused arc, skip
        arc = parentprob.A[a]
        i, j = arc.src, arc.dst
        cost = arc.cost
        M[ia] = max(0, minimum(replace([
            p[ia] - cost,
            omax[j] - omin[i] - cost,
            odmax - (omin[i] + cost + dmin[j]),
            dmax[i] - dmin[j] - cost
        ], NaN => Inf)))
    end

    return M
end

calculate_bigM(prob::EmptyProblem) = zeros(length(tolled_arcs(parent(prob))))
