const AbstractTollsDict = AbstractDict{Int,<:Real}
const TollsDict = Dict{Int,Float64}

zero_tolls(prob::AbstractProblem) = Dict(k => 0.0 for k in tolled_arcs(prob))
inf_tolls(prob::AbstractProblem) = Dict(k => Inf for k in tolled_arcs(prob))