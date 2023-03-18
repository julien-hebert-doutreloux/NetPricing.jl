# Empty problem has no assignment
assign(::Any, ::EmptyProblem; kwargs...) = nothing

# Assign a formulation to a problem type, fallback to standard formulation
# Integrated with improved big-M on PathPreprocessedProblem
function assign(formulation::Type{GeneralFormulation{P,D}}, prob::PathPreprocessedProblem;
    bigM_difference=true, kwargs...) where {P,D}

    form = formulation(prob)
    if bigM_difference
        form = BigMDifference(form, prob)
    end
    return form
end
assign(::Type{GeneralFormulation{P,D}}, prob::AbstractCommodityProblem; kwargs...) where {P,D} = StandardFormulation(prob)
assign(::Type{GeneralFormulation{P,D}}, ::EmptyProblem; kwargs...) where {P,D} = nothing

# Assign formulations based on the number of paths
function assign_breakpoint(probs, formulations, breakpoints; kwargs...)
    forms = Formulation[]
    for prob in probs
        probpaths = paths(prob)
        index = isnothing(probpaths) ? nothing : findfirst(>=(length(probpaths)), breakpoints)
        formtype = isnothing(index) ? StandardFormulation : formulations[index]
        form = assign(formtype, prob; kwargs...)
        isnothing(form) || push!(forms, form)
    end
    return forms
end
