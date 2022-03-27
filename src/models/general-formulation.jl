# Representations
abstract type PrimalRepresentation end
abstract type DualRepresentation end

struct PrimalArc <: PrimalRepresentation
    prob::AbstractCommodityProblem
end
struct PrimalPath <: PrimalRepresentation
    prob::PathPreprocessedProblem
end

struct DualArc <: DualRepresentation
    prob::AbstractCommodityProblem
end
struct DualPath <: DualRepresentation
    prob::PathPreprocessedProblem
end

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
problem(primal::PrimalRepresentation) = primal.prob
problem(dual::DualRepresentation) = dual.prob

# Named formulations
const StandardFormulation = GeneralFormulation{PrimalArc, DualArc}
const PathArcStandardFormulation = GeneralFormulation{PrimalPath, DualArc}
const ValueFunctionFormulation = GeneralFormulation{PrimalArc, DualPath}
const PathValueFunctionFormulation = GeneralFormulation{PrimalPath, DualPath}

const STDFormulation = StandardFormulation
const PASTDFormulation = PathArcStandardFormulation
const VFFormulation = ValueFunctionFormulation
const PVFFormulation = PathValueFunctionFormulation

# General formulation implementation
function Base.append!(model::Model, form::GeneralFormulation, M, N;
    sdtol=1e-10,        # Strong duality tolerance
    kwargs...
    )
    # Primal + linearization
    x, primalobj = formulate_primal(model, form)
    sumtx = linearization(model, form, x, M, N)

    # Dual
    dualobj = formulate_dual(model, form)

    # Strong duality
    @constraint(model, primalobj + sumtx ≤ dualobj + sdtol)

    return sumtx * demand(problem(primal(form)))
end

calculate_bigM(form::GeneralFormulation; kwargs...) = calculate_bigM(problem(primal(form)))
