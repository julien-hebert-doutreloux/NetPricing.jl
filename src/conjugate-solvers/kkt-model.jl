# ConjugateKKTModel: similar to ConjugateLinearModel, but use strong duality to enforce optimality,
# while priority is implemented via objective function.
# This model is able to test for strong bilevel feasibility (see KKTStrongBFTester)
struct ConjugateKKTModel <: AbstractConjugateSolver
    model::Model
    prob::Problem
    num_commodities::Int
    weights::Vector{Float64}
end

# Interface
problem(cmodel::ConjugateKKTModel) = cmodel.prob
num_commodities(cmodel::ConjugateKKTModel) = cmodel.num_commodities
weights(cmodel::ConjugateKKTModel) = cmodel.weights

# Init a conjugate model to test bilevel feasibility between some paths, weights = 1
function ConjugateKKTModel(prob::Problem, num_commodities,
    weights=ones(num_commodities);
    silent=true, threads=nothing, sdtol=default_sdtol())

    model = Model(() -> Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), silent)
    set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)

    nk = num_commodities
    na = length(arcs(prob))
    nv = nodes(prob)
    a1 = tolled_arcs(prob)

    c = cost_vector(prob)
    A = incidence_matrix(prob)
    η = weights

    @variable(model, x[a=1:na, k=1:nk] ≥ 0)
    @variable(model, λ[i=1:nv, k=1:nk])
    @variable(model, t[a=a1] ≥ 0)
    @variable(model, L)                 # The dual objective: L = b'λ - w't
    @variable(model, R)                 # The primal objective (base cost): R = c'x
    @variable(model, s[a=1:na, k=1:nk]) # The slack variables for each dualfeas constraint, is used to test strong bilevel feasibility
    @variable(model, S ≤ 1)             # The global slack variable, connecting individual slack variables

    @objective(model, Max, 0)

    tfull = Array{Any}(zeros(na))
    for a in a1
        tfull[a] = t[a]
    end
    @constraint(model, primalfeas, A * x .== 0)
    @constraint(model, demandlimit[a=a1], η' * @view(x[a,:]) .≤ 0)
    @constraint(model, dualfeas, A' * λ .- tfull .+ s .≤ c)
    @constraint(model, strongdual, L ≥ c' * x * η - sdtol)
    @constraint(model, dualobj, L == 0)
    @constraint(model, slackselect, s .== 0)
    
    return ConjugateKKTModel(model, prob, num_commodities, weights)
end

# Init a conjugate model for the whole problem, weights = commodity demands, odpairs already set
function ConjugateKKTModel(prob::Problem; kwargs...)
    cmodel = ConjugateKKTModel(prob, length(prob.K), [comm.demand for comm in prob.K]; kwargs...)
    set_odpairs(cmodel, [(comm.orig, comm.dest) for comm in prob.K])
    return cmodel
end

# Set OD pairs
function set_odpairs(cmodel::ConjugateKKTModel, odpairs::AbstractVector{Tuple{Int,Int}})
    model, prob, η = cmodel.model, cmodel.prob, cmodel.weights
    b = hcat([sourcesink_vector(prob, o, d) for (o, d) in odpairs]...)
    set_normalized_rhs.(model[:primalfeas], b)
    set_normalized_coefficient.(model[:dualobj], model[:λ], -b .* η')
    return nothing
end

# Set demands (and maximize sum of tolls in priority list)
# If SlackPriority is used, then maximize S instead
struct SlackPriority end

function _set_demands(cmodel::ConjugateKKTModel, demands)
    model = cmodel.model
    for a in tolled_arcs(cmodel.prob)
        w = get(demands, a, 0.)
        set_normalized_coefficient(model[:dualobj], model[:t][a], w)
        set_normalized_rhs(model[:demandlimit][a], w)
    end
    return
end

function set_demands(cmodel::ConjugateKKTModel, demands, ::SlackPriority)
    _set_demands(cmodel, demands)
    @objective(cmodel.model, Max, cmodel.model[:S])
    return
end

function set_demands(cmodel::ConjugateKKTModel, demands, priority)
    _set_demands(cmodel, demands)
    model = cmodel.model
    @objective(model, Max, sum(get(demands, a, 0.) * model[:t][a] for a in priority))
    clear_slack(cmodel)
    return
end

# Optimize
JuMP.optimize!(cmodel::ConjugateKKTModel) = optimize!(cmodel.model)
JuMP.objective_value(cmodel::ConjugateKKTModel) = value(cmodel.model[:L])
tvals(cmodel::ConjugateKKTModel) = value.(cmodel.model[:t]).data

# Functions for strong bilevel feasibility test
function clear_slack(cmodel::ConjugateKKTModel)
    # Set all s to 0
    model = cmodel.model
    set_normalized_coefficient.(model[:slackselect], model[:S], 0)
    return
end

function set_slack(cmodel::ConjugateKKTModel, active_arcs)
    # Connect s to S if the arc is not in active_arcs
    model = cmodel.model
    slackselect, S = model[:slackselect], model[:S]
    na, nk = size(slackselect)
    for k in 1:nk, a in 1:na
        set_normalized_coefficient(slackselect[a, k], S, a ∈ active_arcs[k] ? 0 : -1)
    end
    return
end
