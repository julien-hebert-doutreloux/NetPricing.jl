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
    @constraint(model, sumtx <= unnormalized_objective_term(form) + sdtol)

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

    if isempty(tx)
        sumtx = 0.0
    else
        sumtx = sum(tx)
    end

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

# Envelop only: only adds the constraint, does not make anything binary
struct EnvelopOnly <: CommodityLinearization end

function linearize_commodity_extra(::EnvelopOnly, ::PrimalRepresentation) end



########################## CUSTOM

#### Custom
function custom_linearize!(model::Model, linearization::CommodityLinearization, forms, Ms, N, rtrans, vtrans, ktrans, nv_,na_, c, γc, γA, γt; sdtol=1e-10)
    for (form, M) in zip(forms, Ms)
        custom_linearize_commodity!(model, linearization, form, M, N, rtrans, vtrans, ktrans, nv_,na_, c, γc, γA, γt; sdtol=sdtol)
    end
    return
end

function custom_linearize_commodity!(model::Model, linearization::CommodityLinearization, form::Formulation, M, N, rtrans, vtrans, ktrans, nv_,na_, c, γc, γA, γt; sdtol=1e-10)
    # Linearization
    sumtx = custom_linearize_commodity_primal(model, linearization, primal(form), M, N, rtrans, vtrans, ktrans, nv_,na_, c, γc, γA, γt)
    # Strong duality
    @constraint(model, sumtx <= unnormalized_objective_term(form) + sdtol)
    return sumtx
end

function custom_linearize_commodity_primal(model::Model, linearization::CommodityLinearization, primal::PrimalRepresentation, M, N, rtrans, vtrans, ktrans, nv_,na_, c, γc, γA, γt)

    prob = problem(primal)
    parentprob = parent(prob)

    a1 = tolled_arcs(prob)
    k = index(prob)
    Amap = arcmap(prob)
	Vmap = used_nodes(prob)
	
    x = primal.x[a1]
    tx = @variable(model, [a=a1], lower_bound = 0, base_name="tx[$k]")
    t = remap_t(model, prob)

    if isempty(tx)
        sumtx = 0.0
    else
        sumtx = sum(tx)
    end

    a1dict = Dict(a => i for (i, a) in enumerate(tolled_arcs(parentprob)))
    mapped_a1 = [a1dict[a] for a in Amap[a1]]
    M = @view M[mapped_a1]
    N = @view N[mapped_a1]

    # Linearization
    @constraint(model, tx .≤ M .* x)
    @constraint(model, t .- tx .≥ 0)
    @constraint(model, t .- tx .≤ N .* (1 .- x))
    #########################################################

	
    b = sourcesink_vector(prob)			# Source sink vector 
    nv = length(nodes(prob))			# number of nodes
 	bfull = expand_b(Vmap, nv, b)		# Source sink vector in full dimension
 	println("bfull")
 	γbfull_min, γbfull_avg, γbfull_max   = projection(vtrans, bfull)	# Projected Source sink vector full dimension
 	println("projection")
 	
 	# To verify if the projected problem is feasible
 	if γbfull_min == γbfull_max
 		γbfull = γbfull_min
 	else
 		γbfull = false
 	end
 	
 	if γbfull in ktrans              # it is possible that the problem become infeasible in the transformed space
 				
		k_ = ktrans[projection(vtrans, bfull)]               # retrieve the associated transformed problem
		λ_ = rtrans.λvals[k_]                               # Reduce dual value in transformed space with the corresponding problem
		λ_full = expand_b(rtrans.Vmap[k_], nv_, λ_)         # Full dimension dual value in transformed space
		γ_inv_λ_full = retroprojection(vtrans, value.(λ_full)) # Retroprojection of the dual value into the original space
		
		x_ = rtrans.xvals[k_]
		x_full = expand_b(rtrans.Amap[k_], na_, x_)       # optimal solution path in transformed problem 

		# γ(A)' * λ~ <= γ(c) + γ(t)
		@constraint(model, γA' * λ_full .≤ γc + γt)
		# (c + t)' * x <= b' * γ^-1(λ~)
		@constraint(model, c' * x + sumtx ≤ b' * γ_inv_λ_full)
		# (γ(c) + γ(t))' * x~ <= γ(b)' * λ~
		@constraint(model, (γc + γt)' * x_full ≤ γbfull' * λ_full)
	end
    

    # Run linearization-specific code
    linearize_commodity_extra(linearization, primal)

    return sumtx
end

