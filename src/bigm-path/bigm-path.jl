# Extract big M of each path and overall big M of the arcs
function calculate_bigM_paths(prob::PathPreprocessedProblem, model=inverse_model(parent(prob), 1, silent=true))
    parentprob = parent(prob)
    set_inverse_model_odpairs(model, [index(prob)])    

    paths = prob.original_paths
    np = length(paths)
    a1 = tolled_arcs(parentprob)

    arcdict = srcdst_to_index(parentprob)
    a1set = BitSet(a1)
    a1dict = Dict(a => i for (i, a) in enumerate(a1))

    Mp = zeros(length(a1), np)
    for (p, path) in enumerate(paths)
        tolled = path_tolled_arcs(path, arcdict, a1set)
        for a in tolled
            set_inverse_model_paths(model, [path], Dict(a => 1-1e-5), set_odpairs=false)
            optimize!(model)
            Mp[a1dict[a],p] = value(model[:t][a])
        end
    end

    return Mp, maximum(Mp, dims=2)
end

# Fallback, calculate big M in the old way
function calculate_bigM_paths_fallback(prob::AbstractCommodityProblem, _)
    parentprob = parent(prob)
    graph = build_graph(parentprob)

    a1 = tolled_arcs(parentprob)
    comm = parentprob.K[index(prob)]

    reset_tolls!(graph, parentprob)
    omin = dists_from(graph, comm.orig)
    dmin = dists_to(graph, comm.dest)

    disable_tolls!(graph, parentprob)
    omax = dists_from(graph, comm.orig)
    dmax = dists_to(graph, comm.dest)
    odmax = omax[comm.dest]

    p = [shortest_path(graph, a.src, a.dst)[2] for a in parentprob.A[a1]]

    M = zeros(length(a1))
    for ia in 1:length(a1)
        arc = parentprob.A[a1[ia]]
        i, j = arc.src, arc.dst
        cost = arc.cost
        M[ia] = max(0, minimum(replace([
            p[ia] - cost,
            omax[j] - omin[i] - cost,
            odmax - (omin[i] + cost + dmin[j]),
            dmax[i] - dmin[j] - cost
        ], NaN => Inf)))
    end

    return nothing, M
end

calculate_bigM_paths(prob::EmptyProblem, _) = nothing, zeros(length(tolled_arcs(parent(prob))))
calculate_bigM_paths(prob::AbstractCommodityProblem, model) = calculate_bigM_paths_fallback(prob, model)

function calculate_bigM_paths(prob::PreprocessedProblem, model)
    _, M = calculate_bigM_paths_fallback(prob, model)
    # Set to 0 if arc is not used
    a1 = tolled_arcs(parent(prob))
    a1dict = Dict(a => i for (i, a) in enumerate(a1))
    unused_a1 = setdiff(a1, prob.Amap)
    unused_a1_index = [a1dict[a] for a in unused_a1]
    M[unused_a1_index] .= 0
    return nothing, M
end

# Vector form
function calculate_bigM_paths(probs::Vector{<:AbstractCommodityProblem}; maxpaths=100, silent=true, threads=nothing)
    model = inverse_model(parent(probs[1]), 1; silent=silent, threads=threads)
    probs = [(length(something(paths(prob), [])) > maxpaths ? prob.pprob : prob) for prob in probs]

    nk = length(probs)
    a1 = tolled_arcs(parent(probs[1]))

    Mp = Union{Nothing,Matrix{Float64}}[]
    M = zeros(length(a1), nk)
    for k in 1:nk
        Mpk, Mk = calculate_bigM_paths(probs[k], model)
        push!(Mp, Mpk)
        M[:,k] = Mk
    end

    N = maximum(M, dims=2)[:]

    return Mp, M', N
end