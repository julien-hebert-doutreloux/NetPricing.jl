## Get the big-M parameter for each arc and for each commodity
function calculate_bigM(prob::Problem; graph = build_graph(prob))
    nv = prob.V
    nk = length(prob.K)
    a1 = tolled_arcs(prob)

    omin = zeros(nk, nv)
    omax = zeros(nk, nv)
    dmin = zeros(nk, nv)
    dmax = zeros(nk, nv)

    reset_tolls!(graph, prob)
    for k = 1:nk
        comm = prob.K[k]
        omin[k, :] = dists_from(graph, comm.orig)
        dmin[k, :] = dists_to(graph, comm.dest)
    end

    disable_tolls!(graph, prob)
    for k = 1:nk
        comm = prob.K[k]
        omax[k, :] = dists_from(graph, comm.orig)
        dmax[k, :] = dists_to(graph, comm.dest)
    end

    p = [shortest_path(graph, a.src, a.dst)[2] for a in prob.A[a1]]

    M = zeros(nk, length(a1))
    for k in 1:nk, ia in 1:length(a1)
        arc = prob.A[a1[ia]]
        i, j = arc.src, arc.dst
        cost = arc.cost
        odmax = omax[k, prob.K[k].dest]
        M[k, ia] = max(0, minimum(replace([
            p[ia] - cost,
            omax[k,j] - omin[k,i] - cost,
            odmax - (omin[k,i] + cost + dmin[k,j]),
            dmax[k,i] - dmin[k,j] - cost
        ], NaN => Inf)))
    end

    N = maximum(M, dims=1)'

    return M, N
end

## Big-M calculation using fixed arc solver
function calculate_bigM_fixedarcs(prob::Problem; graph = build_graph(prob), upper=inf_tolls(prob), lower=zero_tolls(prob))
    nv = prob.V
    nk = length(prob.K)
    a1 = tolled_arcs(prob)

    graph2 = copy(graph)
    reset!(graph2, prob)
    set_tolls!(graph2, prob, upper)

    M = zeros(nk, length(a1))
    for k in 1:nk, ia in 1:length(a1)
        a = a1[ia]
        comm = prob.K[k]
        orig, dest = comm.orig, comm.dest

        _, exclude_cost = shortest_path(graph2, orig, dest)
        include_cost = solve_fixedarcs_recur(graph, orig, dest, prob, [a], [], tolls=lower)

        M[k, ia] = max(0, exclude_cost - include_cost + lower[a])
    end

    N = maximum(M, dims=1)'

    return M, N
end

function calculate_bigM_iterative(prob::Problem, iterations=10; graph = build_graph(prob), upper=inf_tolls(prob), lower=zero_tolls(prob))
    a1 = tolled_arcs(prob)

    local M, N

    for t in 1:iterations
        M, N = calculate_bigM_fixedarcs(prob, graph=graph, upper=upper, lower=lower)
        upper = Dict(zip(a1, N))
    end

    return M, N
end

initial_tmax(prob::Problem) = Dict(zip(tolled_arcs(prob), calculate_bigM(prob)[2]))