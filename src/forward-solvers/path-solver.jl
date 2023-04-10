mutable struct CommodityForwardPathSolver <: CommodityForwardSolver
    prob::PathPreprocessedProblem
    ts::Vector{Float64}
    optimal_a1::BitSet
    optimal_cost::Float64
    parent_a1::Vector{Int}
    pc::Vector{Float64}
    pa1::Vector{BitSet}
end

problem(fmodel::CommodityForwardPathSolver) = fmodel.prob

function CommodityForwardPathSolver(prob::PathPreprocessedProblem)
    na = length(arcs(parent(prob)))
    parent_a1 = tolled_arcs(parent(prob))
    ts = zeros(na)

    a1set = BitSet(tolled_arcs(prob))
    arcdict = srcdst_to_index(prob)
    arccost = srcdst_to_cost(prob)
    Amap = used_arcs(prob)

    ps = prob.paths
    pc = get_path_cost.(ps, Ref(arccost))
    pa1 = [BitSet(Amap[path_tolled_arcs(p, arcdict, a1set)]) for p in ps]

    return CommodityForwardPathSolver(prob, ts, BitSet(), 0., parent_a1, pc, pa1)
end

function set_tolls(fmodel::CommodityForwardPathSolver, tolls)
    ts = fmodel.ts
    for a in fmodel.parent_a1
        ts[a] = tolls[a]
    end
    return nothing
end

function JuMP.optimize!(fmodel::CommodityForwardPathSolver)
    ts = fmodel.ts
    best_obj, best_pa1 = Inf, BitSet()
    for (cp, a1p) in zip(fmodel.pc, fmodel.pa1)
        (cp >= best_obj) && break
        cost = cp + sum(ts[a] for a in a1p; init = 0.)
        if cost < best_obj
            best_obj, best_pa1 = cost, a1p
        end
    end
    fmodel.optimal_cost = best_obj
    fmodel.optimal_a1 = best_pa1
end

JuMP.objective_value(fmodel::CommodityForwardPathSolver) = fmodel.optimal_cost
optimal_tolled_set(fmodel::CommodityForwardPathSolver) = fmodel.optimal_a1
