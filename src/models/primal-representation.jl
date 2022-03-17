# Primal arc
function primal_arc(model::Model, prob::AbstractCommodityProblem)
    na = length(arcs(prob))
    a1 = tolled_arcs(prob)
    k = index(prob)
    
    c = cost_vector(prob)
    A = incidence_matrix(prob)
    b = sourcesink_vector(prob)

    x = @variable(model, [a=1:na], lower_bound = 0, upper_bound = 1, binary = a in a1, base_name="x[$k]")

    primalfeas = @constraint(model, A * x .== b)
    primalobj = c' * x

    return x, primalfeas, primalobj, nothing, nothing
end

# Primal path
primal_path(model::Model, prob::AbstractCommodityProblem) = primal_arc(model, prob)

function primal_path(model::Model, prob::PathPreprocessedProblem)
    all_paths = prob.paths
    np = length(all_paths)
    a1 = tolled_arcs(prob)
    k = index(prob)

    z = @variable(model, [1:np], lower_bound = 0, upper_bound = 1, binary = true, base_name="z[$k]")
    x = @variable(model, [a=a1], lower_bound = 0, upper_bound = 1, base_name="x[$k]")

    primalfeas = @constraint(model, sum(z) == 1)
    conversion = @constraint(model, [a=a1], 0 == x[a])

    for (p, path) in enumerate(all_paths)
        for a in path_arcs(path, prob)
            (a in a1) || continue
            set_normalized_coefficient(conversion[a], z[p], 1)
        end
    end

    arccosts = srcdst_to_cost(prob)
    pathcosts = [get_path_cost(path, arccosts) for path in all_paths]
    primalobj = pathcosts' * z

    return x, primalfeas, primalobj, z, conversion
end

# Linearization
function linearization(model::Model, prob::AbstractCommodityProblem, x, M, N)
    a1 = tolled_arcs(prob)
    k = index(prob)
    Amap = arcmap(prob)

    tx = @variable(model, [a=a1], lower_bound = 0, base_name="tx[$k]")
    t = JuMP.Containers.DenseAxisArray(model[:t][Amap[a1]].data, a1)

    sumtx = sum(tx)

    a1dict = Dict(a => i for (i, a) in enumerate(tolled_arcs(parent(prob))))
    mapped_a1 = [a1dict[a] for a in Amap[a1]]
    M = @view M[mapped_a1]
    N = @view N[mapped_a1]

    bilinear1 = @constraint(model, tx .≥ 0)
    bilinear2 = @constraint(model, tx .≤ M .* x[a1])
    bilinear3 = @constraint(model, t .- tx .≥ 0)
    bilinear4 = @constraint(model, t .- tx .≤ N .* (1 .- x[a1]))

    return tx, sumtx, bilinear1, bilinear2, bilinear3, bilinear4
end