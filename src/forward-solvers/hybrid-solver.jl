struct ForwardHybridSolver <: AbstractForwardSolver
    solvers::Vector{CommodityForwardSolver}
    weights::Vector{Float64}
    a1::Vector{Int}
end

problem(fmodel::ForwardHybridSolver) = fmodel.solvers |> first |> problem |> parent

function ForwardHybridSolver(probs::Vector{<:AbstractCommodityProblem}, weights=demand.(probs); maxpaths = 1000)
    solvers = _assign_forward_solver.(probs; maxpaths)
    return ForwardHybridSolver(solvers, weights, tolled_arcs(parent(first(probs))))
end

_assign_forward_solver(prob::EmptyProblem; kwargs...) = CommodityForwardEmptySolver(prob)
_assign_forward_solver(prob::AbstractCommodityProblem; kwargs...) = CommodityForwardGraphSolver(prob)
_assign_forward_solver(prob::PathPreprocessedProblem; maxpaths) =
    length(paths(prob)) <= maxpaths ? CommodityForwardPathSolver(prob) : CommodityForwardGraphSolver(prob)

function set_tolls(fmodel::ForwardHybridSolver, tolls, priority; discount=default_discount())
    toll_dict = Dict(a => get(tolls, a, 0.) * (a ∈ priority ? discount : 1.) for a in fmodel.a1)
    for solver in fmodel.solvers
        set_tolls(solver, toll_dict)
    end
    return nothing
end

JuMP.optimize!(fmodel::ForwardHybridSolver) = optimize!.(fmodel.solvers)

JuMP.objective_value(fmodel::ForwardHybridSolver) = sum(objective_value(solver) * w for (solver, w) in zip(fmodel.solvers, fmodel.weights))

function wvals(fmodel::ForwardHybridSolver)
    na = fmodel |> problem |> arcs |> length
    ws = zeros(na)
    for (solver, w) in zip(fmodel.solvers, fmodel.weights)
        pa1 = optimal_tolled_set(solver)
        for a in pa1
            ws[a] += w
        end
    end
    return ws[fmodel.a1]
end
