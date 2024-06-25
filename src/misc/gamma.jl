function expand_b(Vmap, nv, b)
    # b is the vector to expand in the sparse space of dimension nv
    # Vmap is the index mapping to place element of b in the vector of dimension nv
    bfull = Vector{AffExpr}(zeros(nv))
    bfull[Vmap] .= b
    return bfull
end

function projection(transformation::Dict, Nn)
    NS = [[] for _ in 1:maximum(values(transformation))]
    for i in keys(transformation)
    	if typeof(i)!=typeof(1)
    		j = parse(Int, i)
    	else
    		j=i
    	end
        append!(
                NS[transformation[i]], Nn[j]#N[parse(Int, i)]
        )
    end
    
    NT_min = [value(minimum(vec)) for vec in NS]
    NT_avg = [value(mean(vec)) for vec in NS]
    NT_max = [value(maximum(vec)) for vec in NS]
    return NT_min, NT_avg, NT_max
end
function retroprojection(transformation::Dict, NT)
    # transformation is a dict of the for "i":j
    Nn = zeros(length(transformation))
    for (k, v) in transformation
    	if typeof(k)!=typeof(1)
    		j = parse(Int, k)
    	else
    		j=k
    	end
        Nn[j] = NT[v]
    end
    
    return Nn
end


