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

# Query
value_x(model) = value.(model[:x])
value_t(model) = value.(model[:t])

# Fix var
fix_var(var::VariableRef) = fix(var, 0, force=true)