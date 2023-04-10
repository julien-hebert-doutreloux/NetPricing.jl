mutable struct CommodityForwardGraphSolver <: CommodityForwardSolver
    graph::AbstractGraph
    prob::AbstractCommodityProblem
    optimal_a1::BitSet
    optimal_cost::Float64
    a1set::BitSet
    arcdict::Dict
end

problem(fmodel::CommodityForwardGraphSolver) = fmodel.prob

function CommodityForwardGraphSolver(prob::AbstractCommodityProblem)
    graph = build_graph(prob)
    a1set = BitSet(tolled_arcs(prob))
    arcdict = srcdst_to_index(prob)
    return CommodityForwardGraphSolver(graph, prob, BitSet(), 0., a1set, arcdict)
end

function set_tolls(fmodel::CommodityForwardGraphSolver, tolls)
    set_tolls!(fmodel.graph, fmodel.prob, tolls)
    return nothing
end

function JuMP.optimize!(fmodel::CommodityForwardGraphSolver)
    path, fmodel.optimal_cost = shortest_path(fmodel.graph, orig(fmodel.prob), dest(fmodel.prob))
    pa1 = path_tolled_arcs(path, fmodel.arcdict, fmodel.a1set)
    fmodel.optimal_a1 = BitSet(used_arcs(fmodel.prob)[pa1])
end

JuMP.objective_value(fmodel::CommodityForwardGraphSolver) = fmodel.optimal_cost
optimal_tolled_set(fmodel::CommodityForwardGraphSolver) = fmodel.optimal_a1
