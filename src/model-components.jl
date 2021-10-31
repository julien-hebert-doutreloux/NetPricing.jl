# Cost vector
cost_vector(prob::Problem) = [arc.cost for arc in prob.A]

# Incidence matrix
function incidence_matrix(prob::Problem)
    mat = spzeros(Int, prob.V, length(prob.A))
    for (i, arc) in enumerate(prob.A)
        mat[arc.src, i] = 1
        mat[arc.dst, i] = -1
    end
    return mat
end

# Source-sink vector
function sourcesink_vector(prob::Problem, o, d)
    b = spzeros(Int, prob.V)
    b[o] = 1
    b[d] = -1
    return b
end

sourcesink_vector(prob::Problem, k) = sourcesink_vector(prob, prob.K[k].orig, prob.K[k].dest)

# Demand vector
demand_vector(prob::Problem) = [comm.demand for comm in prob.K]

# Query
value_x(model) = value.(model[:x])
value_t(model) = value.(model[:t])