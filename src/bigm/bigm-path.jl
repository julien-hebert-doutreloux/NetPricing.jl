# Extract big M of each path and overall big M of the arcs
# ConjugateLinearModel is faster than ConjugateDynamicLinearModel if there is only 1 commodity
function calculate_bigM_paths(prob::PathPreprocessedProblem; threads=nothing, conjugate_solver=ConjugateLinearModel)
    parentprob = parent(prob)
    # 
    cmodel = conjugate_solver(parent(prob), 1, threads=threads)
    set_odpairs(cmodel, [index(prob)])

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
            set_paths(cmodel, [path], Dict(a => 1-1e-5), set_odpairs=false)
            optimize!(cmodel)
            Mp[a1dict[a],p] = value(cmodel.model[:t][a])
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
primal(form::BigMPath) = primal(form.form)
dual(form::BigMPath) = dual(form.form)
objective_term(form::BigMPath) = objective_term(form.form)

Base.append!(model::Model, form::BigMPath; kwargs...) = append!(model, form.form; kwargs...)

calculate_bigM(form::BigMPath; threads=nothing, conjugate_solver=ConjugateLinearModel, kwargs...) =
    calculate_bigM_paths(form.prob, threads=threads, conjugate_solver=conjugate_solver)

Base.show(io::IO, form::BigMPath) = print(io, form.form, " + Big-M from paths")
