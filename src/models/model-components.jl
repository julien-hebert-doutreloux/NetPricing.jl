# Cost vector
cost_vector(prob::AbstractProblem) = [arc.cost for arc in arcs(prob)]

# Incidence matrix
function incidence_matrix(prob::AbstractProblem)
    mat = spzeros(Int, nodes(prob), length(arcs(prob)))
    for (i, arc) in enumerate(arcs(prob))
        mat[arc.src, i] = 1
        mat[arc.dst, i] = -1
    end
    return mat
end

# Source-sink vector
function sourcesink_vector(prob::AbstractProblem, o, d)
    b = zeros(Int, nodes(prob))
    b[o] = 1
    b[d] = -1
    return b
end

sourcesink_vector(prob::Problem, k) = sourcesink_vector(prob, prob.K[k].orig, prob.K[k].dest)
sourcesink_vector(prob::AbstractCommodityProblem) = sourcesink_vector(prob, orig(prob), dest(prob))
sourcesink_vector(::EmptyProblem) = Int[]

# Demand vector
demand_vector(prob::Problem) = [comm.demand for comm in prob.K]

# Path-arc incident
function path_arc_incident_matrix(prob::PathPreprocessedProblem)
    mat = spzeros(Int, length(arcs(prob)), length(prob.paths))
    arcdict = srcdst_to_index(prob)
    for (p, path) in enumerate(prob.paths)
        for a in path_arcs(path, arcdict)
            mat[a, p] = 1
        end
    end
    return mat
end

# Utils
function remap_t(model, prob::AbstractCommodityProblem)
    Amap = arcmap(prob)
    a1 = tolled_arcs(prob)
    return DenseAxisArray(model[:t][Amap[a1]].data, a1)
end

function expand_t(t, prob::AbstractCommodityProblem)
    na = length(arcs(prob))
    a1 = tolled_arcs(prob)
    tfull = Vector{AffExpr}(zeros(na))
    tfull[a1] .= t.data
    return tfull
end

# Query
value_x(model) = value.(model[:x])
value_t(model) = value.(model[:t])

# Fix var
fix_var(var::VariableRef) = fix(var, 0, force=true)