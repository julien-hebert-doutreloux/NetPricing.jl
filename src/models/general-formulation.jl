# General Formulation with selectable primal and dual
struct GeneralFormulation{P<:PrimalRepresentation, D<:DualRepresentation} <: Formulation
    primal::P
    dual::D
end

GeneralFormulation{P,D}(prob) where {P,D} = GeneralFormulation(P(prob), D(prob))

# Type queries
primal_type(::Type{GeneralFormulation{P,D}}) where {P,D} = P
dual_type(::Type{GeneralFormulation{P,D}}) where {P,D} = D

primal_type(form::GeneralFormulation) = primal_type(typeof(form))
dual_type(form::GeneralFormulation) = dual_type(typeof(form))

# Accessors
primal(form::GeneralFormulation) = form.primal
dual(form::GeneralFormulation) = form.dual

problem(form::GeneralFormulation) = parent(problem(primal(form)))

# General formulation implementation
function Base.append!(model::Model, form::GeneralFormulation; kwargs...)
    formulate_primal!(model, primal(form))
    formulate_dual!(model, dual(form))
    return
end

calculate_bigM(form::GeneralFormulation; kwargs...) = calculate_bigM(problem(primal(form)))
objective_term(form::GeneralFormulation) = demand(problem(primal(form))) * (dualobj(form) - primalobj(form))

# Pretty print
Base.show(io::IO, form::GeneralFormulation{P,D}) where {P,D} = print(io, "GeneralFormulation{$P, $D}(", problem(primal(form)), ")")