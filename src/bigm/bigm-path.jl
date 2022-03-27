# Extract big M of each path and overall big M of the arcs
function calculate_bigM_paths(prob::PathPreprocessedProblem; threads=nothing)
    parentprob = parent(prob)
    model = inverse_model(parent(prob), 1, silent=true, threads=threads)
    set_inverse_model_odpairs(model, [index(prob)])

    paths = prob.original_paths
    np = length(paths)
    a1 = tolled_arcs(parentprob)

    arcdict = srcdst_to_index(parentprob)
    a1set = BitSet(a1)
    a1dict = Dict(a => i for (i, a) in enumerate(a1))

    Mp = spzeros(length(a1), np)
    for (p, path) in enumerate(paths)
        tolled = path_tolled_arcs(path, arcdict, a1set)
        for a in tolled
            set_inverse_model_paths(model, [path], Dict(a => 1-1e-5), set_odpairs=false)
            optimize!(model)
            Mp[a1dict[a],p] = value(model[:t][a])
        end
    end

    return maximum(Mp, dims=2)[:]
end

# Wrapper to mark that this formulation use calculate_bigM_paths
struct BigMPath{T<:Formulation} <: Formulation
    form::T
    prob::PathPreprocessedProblem
end

problem(form::BigMPath) = problem(form.form)

Base.append!(model::Model, form::BigMPath, args...; kwargs...) = append!(model, form.form, args...; kwargs...)
calculate_bigM(form::BigMPath; threads=nothing, kwargs...) = calculate_bigM_paths(form.prob, threads=threads)

Base.show(io::IO, form::BigMPath) = print(io, form.form, " + Big-M from paths")
