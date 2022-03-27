function general_model(form_type::Type{GeneralFormulation{P,D}}, probs; kwargs...) where {P,D}
    forms = assign.(form_type, probs)
    return formulate(forms; kwargs...)
end

standard_model(probs; kwargs...) = general_model(StandardFormulation, probs; kwargs...)
const std_model = standard_model

path_arc_standard_model(probs; kwargs...) = general_model(PathArcStandardFormulation, probs; kwargs...)
const pastd_model = path_arc_standard_model

value_function_model(probs; kwargs...) = general_model(ValueFunctionFormulation, probs; kwargs...)
const vf_model = value_function_model

path_value_function_model(probs; kwargs...) = general_model(PathValueFunctionFormulation, probs; kwargs...)
const pvf_model = path_value_function_model
