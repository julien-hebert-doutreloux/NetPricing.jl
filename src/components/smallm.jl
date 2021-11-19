## Big-M calculation using fixed arc solver
function calculate_smallM_fixedarcs(prob::Problem; graph = build_graph(prob), upper=inf_tolls(prob), lower=zero_tolls(prob))
    nv = prob.V
    nk = length(prob.K)
    a1 = tolled_arcs(prob)

    M = zeros(nk, length(a1))
    for k in 1:nk, ia in 1:length(a1)
        a = a1[ia]
        comm = prob.K[k]
        orig, dest = comm.orig, comm.dest

        exclude_cost = solve_fixedarcs_hungarian(graph, orig, dest, prob, [], [a], tolls=lower)
        include_cost = solve_fixedarcs_hungarian(graph, orig, dest, prob, [a], [], tolls=upper, relaxed=true)

        M[k, ia] = max(0, exclude_cost - include_cost + lower[a])
    end

    N = minimum(M, dims=1)'

    return M, N
end
