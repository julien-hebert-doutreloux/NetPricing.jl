# Named formulations
const StandardFormulation = GeneralFormulation{PrimalArc, DualArc}
const PathArcStandardFormulation = GeneralFormulation{PrimalPath, DualArc}
const ValueFunctionFormulation = GeneralFormulation{PrimalArc, DualPath}
const PathValueFunctionFormulation = GeneralFormulation{PrimalPath, DualPath}

const STDFormulation = StandardFormulation
const PASTDFormulation = PathArcStandardFormulation
const VFFormulation = ValueFunctionFormulation
const PVFFormulation = PathValueFunctionFormulation

# Named models
function general_model(form_type::Type{GeneralFormulation{P,D}}, probs;
    bigM_difference=true,
    linearization=ArcLinearization(),
    kwargs...) where {P,D}

    forms = convert.(Formulation, filter!(!isnothing, assign.(form_type, probs; bigM_difference)))
    return formulate!(forms, linearization; kwargs...)
end

standard_model(probs; kwargs...) = general_model(StandardFormulation, probs; kwargs...)
const std_model = standard_model

path_arc_standard_model(probs; kwargs...) = general_model(PathArcStandardFormulation, probs; kwargs...)
const pastd_model = path_arc_standard_model

value_function_model(probs; kwargs...) = general_model(ValueFunctionFormulation, probs; kwargs...)
const vf_model = value_function_model

path_value_function_model(probs; kwargs...) = general_model(PathValueFunctionFormulation, probs; kwargs...)
const pvf_model = path_value_function_model



# Custom
function custom_model(form_type::Type{GeneralFormulation{P,D}}, probs, M, N; option=0,
    linearization=ArcLinearization(),
    kwargs...) where {P,D}

    forms = convert.(Formulation, filter!(!isnothing, assign.(form_type, probs)))
    return custom_formulate!(forms, linearization, M, N, option; kwargs...)
end

custom_standard_model(probs, M, N;option, kwargs...) = custom_model(StandardFormulation, probs, M, N;option, kwargs...)
const cst_model = custom_standard_model
