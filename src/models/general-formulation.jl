# General Formulation with selectable primal and dual
struct GeneralFormulation{P<:PrimalRepresentation, D<:DualRepresentation} <: Formulation
    primal::P
    dual::D
end

GeneralFormulation{P,D}(prob; binary_x=false) where {P,D} = GeneralFormulation(P(prob; binary_x=binary_x), D(prob))

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
function Base.append!(model::Model, form::GeneralFormulation, M, N;
    sdtol=1e-10,        # Strong duality tolerance
    kwargs...
    )
    # Primal + linearization
    x, primalobj = formulate_primal!(model, primal(form))
    sumtx = linearization(model, problem(primal(form)), x, M, N)

    # Dual
    dualobj = formulate_dual!(model, dual(form))

    # Strong duality
    @constraint(model, primalobj + sumtx ≤ dualobj + sdtol)

    return sumtx * demand(problem(primal(form)))
end

calculate_bigM(form::GeneralFormulation; kwargs...) = calculate_bigM(problem(primal(form)))

# Pretty print
Base.show(io::IO, form::GeneralFormulation{P,D}) where {P,D} = print(io, "GeneralFormulation{$P, $D}(", problem(primal(form)), ")")