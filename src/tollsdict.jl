const AbstractTollsDict = AbstractDict{Int,<:Real}
const TollsDict = Dict{Int,Float64}

zero_tolls(prob::Problem) = Dict(k => 0.0 for k in tolled_arcs(prob))
inf_tolls(prob::Problem) = Dict(k => Inf for k in tolled_arcs(prob))