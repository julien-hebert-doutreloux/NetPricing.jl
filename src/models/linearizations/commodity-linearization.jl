# Commodity linearization: linearize by commodity
abstract type CommodityLinearization <: AbstractLinearization end

function linearize!(model::Model, linearization::CommodityLinearization, forms, Ms, N; sdtol=1e-10)
    for (form, M) in zip(forms, Ms)
        linearize_commodity!(model, linearization, form, M, N, sdtol=sdtol)
    end
    return
end

function linearize_commodity!(model::Model, linearization::CommodityLinearization, form::Formulation, M, N; sdtol=1e-10)
    # Linearization
    sumtx = linearize_commodity_primal(model, linearization, primal(form), M, N)

    # Strong duality
    @constraint(model, primalobj(form) + sumtx <= dualobj(form) + sdtol)

    return sumtx
end

function linearize_commodity_primal(model::Model, linearization::CommodityLinearization, primal::PrimalRepresentation, M, N)
    prob = problem(primal)
    parentprob = parent(prob)

    a1 = tolled_arcs(prob)
    k = index(prob)
    Amap = arcmap(prob)

    x = primal.x[a1]
    tx = @variable(model, [a=a1], lower_bound = 0, base_name="tx[$k]")
    t = remap_t(model, prob)

    sumtx = sum(tx)

    a1dict = Dict(a => i for (i, a) in enumerate(tolled_arcs(parentprob)))
    mapped_a1 = [a1dict[a] for a in Amap[a1]]
    M = @view M[mapped_a1]
    N = @view N[mapped_a1]

    # Linearization
    @constraint(model, tx .≤ M .* x)
    @constraint(model, t .- tx .≥ 0)
    @constraint(model, t .- tx .≤ N .* (1 .- x))

    # Run linearization-specific code
    linearize_commodity_extra(linearization, primal)

    return sumtx
end

# Arc linearization: make x binary
struct ArcLinearization <: CommodityLinearization end

function linearize_commodity_extra(::ArcLinearization, primal::PrimalRepresentation)
    prob = problem(primal)
    a1 = tolled_arcs(prob)
    x = primal.x
    set_binary.(x[a1])
end

# Path linearization: make z binary if available
struct PathLinearization <: CommodityLinearization end

function linearize_commodity_extra(::PathLinearization, primal::PrimalPath)
    set_binary.(primal.z)
end

linearize_commodity_extra(::PathLinearization, primal::PrimalArc) =
    linearize_commodity_extra(ArcLinearization(), primal)
