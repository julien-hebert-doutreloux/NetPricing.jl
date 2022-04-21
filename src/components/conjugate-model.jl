struct ConjugateModel
    model::Model
    prob::Problem
    num_commodities::Int
    weights::Vector{Float64}
end

# Init a conjugate model to test bilevel feasibility between some paths, weights = 1
function ConjugateModel(prob::Problem, num_commodities,
    weights=ones(num_commodities);
    silent=false, threads=nothing)

    model = Model(() -> Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), silent)
    set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)

    nk = num_commodities
    na = length(arcs(prob))
    nv = nodes(prob)
    a1 = tolled_arcs(prob)

    c = cost_vector(prob)
    A = incidence_matrix(prob)

    @variable(model, λ[i=1:nv, k=1:nk])
    @variable(model, t[a=a1] ≥ 0)

    @objective(model, Max, 0)

    tfull = Array{Any}(zeros(na))
    for a in a1
        tfull[a] = t[a]
    end
    @constraint(model, dualfeas, A' * λ .- tfull .≤ c)
    
    return ConjugateModel(model, prob, num_commodities, weights)
end

# Init a conjugate model for the whole problem, weights = commodity demands, odpairs already set
function ConjugateModel(prob::Problem; kwargs...)
    cmodel = ConjugateModel(prob, length(prob.K), [comm.demand for comm in prob.K]; kwargs...)
    set_odpairs(cmodel, [(comm.orig, comm.dest) for comm in prob.K])
    return cmodel
end

# Wrapper
JuMP.optimize!(cmodel::ConjugateModel) = optimize!(cmodel.model)
JuMP.objective_value(cmodel::ConjugateModel) = objective_value(cmodel.model)
tvals(cmodel::ConjugateModel) = value.(cmodel.model[:t])

# Set OD pairs
function set_odpairs(cmodel::ConjugateModel, odpairs::AbstractVector{Tuple{Int,Int}})
    model, prob, η = cmodel.model, cmodel.prob, cmodel.weights
    b = hcat([sourcesink_vector(prob, o, d) for (o, d) in odpairs]...) .* η'
    set_objective_coefficient.(model, model[:λ], b)
    return nothing
end

set_odpairs(cmodel::ConjugateModel, commodities::AbstractVector{Commodity}) =
    set_odpairs(cmodel, [(comm.orig, comm.dest) for comm in commodities])

set_odpairs(cmodel::ConjugateModel, ks::AbstractVector{Int}) =
    set_odpairs(cmodel, cmodel.prob.K[ks])

# Set demands
function set_demands(cmodel::ConjugateModel, demands, discount::AbstractDict{Int,<:Real})
    model = cmodel.model
    for (a, w) in demands
        set_objective_coefficient(model, model[:t][a], -w * get(discount, a, 1))
    end
    return nothing
end

function set_demands(cmodel::ConjugateModel, demands, discount::Real = 1-1e-5)
    a1 = tolled_arcs(cmodel.prob)
    set_demands(cmodel, demands, Dict(a1 .=> discount))
end

# Set paths
function set_paths(cmodel::ConjugateModel, paths, discount = 1-1e-5; set_odpairs=true)
    prob = cmodel.prob
    nk, η = cmodel.num_commodities, cmodel.weights
    a1 = tolled_arcs(prob)

    # Set od pairs if required
    if set_odpairs
        NetPricing.set_odpairs(cmodel, [(path[1], path[end]) for path in paths])
    end

    # Extract the set of tolled arcs for each path
    arcdict = srcdst_to_index(prob)
    a1set = BitSet(a1)
    tolled_sets = [BitSet(path_tolled_arcs(path, arcdict, a1set)) for path in paths]

    # Set the demands
    demands = Dict(a => sum(a in tolled_sets[k] ? η[k] : 0 for k in 1:nk) for a in a1)
    set_demands(cmodel, demands, discount)
end

# Bilevel feasible test
function is_bilevel_feasible(cmodel::ConjugateModel, paths; kwargs...)
    set_paths(cmodel, paths, 1; kwargs...)  # The value of t is not important, hence no discount
    optimize!(cmodel)
    objval = objective_value(cmodel)
    arccosts = srcdst_to_cost(cmodel.prob)
    expecting = sum(get_path_cost(path, arccosts) for path in paths)
    return objval >= expecting - 1e-6
end

function is_bilevel_feasible(prob::Problem, paths; threads=nothing)
    cmodel = ConjugateModel(prob, length(paths), silent=true, threads=threads)
    return is_bilevel_feasible(cmodel, paths)
end