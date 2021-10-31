## Get the big-M parameter for each arc and for each commodity
function calculate_bigM(prob::Problem; graph = build_graph(prob))
    nv = prob.V
    nk = length(prob.K)
    a1 = tolled_arcs(prob)

    omin = zeros(nk, nv)
    omax = zeros(nk, nv)
    dmin = zeros(nk, nv)
    dmax = zeros(nk, nv)

    set_tolls!(graph, prob, zero_tolls(prob))
    for k = 1:nk, v = 1:nv
        comm = prob.K[k]
        omin[k, v] = shortest_path(graph, comm.orig, v)[2]
        dmin[k, v] = shortest_path(graph, v, comm.dest)[2]
    end

    set_tolls!(graph, prob, inf_tolls(prob))
    for k = 1:nk, v = 1:nv
        comm = prob.K[k]
        omax[k, v] = shortest_path(graph, comm.orig, v)[2]
        dmax[k, v] = shortest_path(graph, v, comm.dest)[2]
    end

    p = [shortest_path(graph, a.src, a.dst)[2] for a in prob.A[a1]]

    M = zeros(nk, length(a1))
    for k in 1:nk, ia in 1:length(a1)
        arc = prob.A[a1[ia]]
        i, j = arc.src, arc.dst
        cost = arc.cost
        odmax = omax[k, prob.K[k].dest]
        M[k, ia] = max(0, min(
            p[ia] - cost,
            omax[k,j] - omin[k,i] - cost,
            odmax - (omin[k,i] + cost + dmin[k,j]),
            dmax[k,i] - dmin[k,j] - cost))
    end

    N = maximum(M, dims=1)'

    return M, N
end

initial_tmax(prob::Problem) = Dict(zip(tolled_arcs(prob), calculate_bigM(prob)[2]))