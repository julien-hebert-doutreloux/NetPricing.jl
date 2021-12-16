function preprocess_path(prob::Problem, k, paths)
    arcindices = srcdst_to_index(prob)

    used_arcs = BitSet()
    used_nodes = BitSet()
    for path in paths
        for i in 1:(length(path)-1)
            push!(used_arcs, arcindices[(path[i], path[i+1])])
            push!(used_nodes, path[i])
        end
        push!(used_nodes, path[end])
    end

    Vmap = collect(used_nodes)
    Amap = collect(used_arcs)

    Vrevmap = revmap(Vmap, prob.V)
    Arevmap = revmap(Amap, length(prob.A))

    # Translate the arcs
    A = [ProblemArc(Vrevmap[arc.src], Vrevmap[arc.dst], arc.cost, arc.toll) for arc in @view prob.A[Amap]]

    # Translate the paths
    mapped_paths = [Vrevmap[path] for path in paths]

    return PathPreprocessedProblem(
        PreprocessedProblem(prob, length(Vmap), A, Vmap, Amap, Vrevmap, Arevmap, k,
        Vrevmap[prob.K[k].orig],
        Vrevmap[prob.K[k].dest]),
        mapped_paths,
        paths
    )
end
