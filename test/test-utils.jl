islessapprox(a, b; kwargs...) = (a <= b) || isapprox(a, b; kwargs...)
≲ = islessapprox
