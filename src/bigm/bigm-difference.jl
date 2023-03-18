# Extract big M using the difference of semi-conjugate functions
function calculate_bigM_difference(prob::PathPreprocessedProblem)
    parentprob = parent(prob)
    Vmap = used_nodes(prob)
    o, d = Vmap[orig(prob)], Vmap[dest(prob)]

    arcdict = srcdst_to_index(parentprob)
    arccosts = srcdst_to_cost(parentprob)
    a1 = tolled_arcs(parentprob)
    a1set = BitSet(a1)
    a1dict = Dict(a => i for (i, a) in enumerate(a1))

    node_paths = prob.original_paths
    tolled_arc_paths = path_tolled_arcs.(node_paths, Ref(arcdict), Ref(a1set))
    paths = collect(zip(node_paths, tolled_arc_paths))

    graph = build_graph(parentprob)
    init_ws = copy(graph.weights)

    M = zeros(length(a1))

    for a in a1
        paths_a = filter(p -> a in p[2], paths)
        isempty(paths_a) && continue
        for (node_p, arc_p) in paths_a
            disabled = collect(setdiff(a1set, arc_p) ∪ a)
            disable_arcs!(graph, parentprob, disabled)
            _, cost0 = shortest_path(graph, o, d)
            restore!(graph, init_ws)
            cost1 = get_path_cost(node_p, arccosts)
            M[a1dict[a]] = max(M[a1dict[a]], cost0 - cost1)
        end
    end

    return M
end

# Wrapper to mark that this formulation use calculate_bigM_paths
struct BigMDifference{T<:Formulation} <: Formulation
    form::T
    prob::PathPreprocessedProblem
end

problem(form::BigMDifference) = problem(form.form)
primal(form::BigMDifference) = primal(form.form)
dual(form::BigMDifference) = dual(form.form)
objective_term(form::BigMDifference) = objective_term(form.form)
unnormalized_objective_term(form::BigMDifference) = unnormalized_objective_term(form.form)

Base.append!(model::Model, form::BigMDifference; kwargs...) = append!(model, form.form; kwargs...)

calculate_bigM(form::BigMDifference; kwargs...) = calculate_bigM_difference(form.prob)

Base.show(io::IO, form::BigMDifference) = print(io, form.form, " + Big-M from cost difference")
