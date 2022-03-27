# Empty problem has no assignment
assign(::Any, ::EmptyProblem) = nothing

# Assign a formulation to a problem type, fallback to standard formulation
assign(formulation::Type{GeneralFormulation{P,D}}, prob::PathPreprocessedProblem) where {P,D} = formulation(prob)
assign(::Type{GeneralFormulation{P,D}}, prob::AbstractCommodityProblem) where {P,D} = StandardFormulation(prob)

# Assign formulations based on the number of paths
function assign_breakpoint(probs, formulations, breakpoints)
    forms = Formulation[]
    for prob in probs
        probpaths = paths(prob)
        index = isnothing(probpaths) ? nothing : findfirst(>=(length(probpaths)), breakpoints)
        form = isnothing(index) ? StandardFormulation : formulations[index]
        push!(forms, assign(form, prob))
    end
    # Filter out all empty formulation
    filter!(!isnothing, forms)
    return forms
end
