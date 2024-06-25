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
function custom_linearize!(model::Model, linearization::CommodityLinearization, forms, Ms, N, rtrans, vtrans, ktrans, nv, nv_, na_, c, γc, γA, γt, γa1; sdtol=1e-10)
    for (form, M) in zip(forms, Ms)
        custom_linearize_commodity!(model, linearization, form, M, N, rtrans, vtrans, ktrans, nv, nv_, na_, c, γc, γA, γt, γa1; sdtol=sdtol)
    end
    return
end

function custom_linearize_commodity!(model::Model, linearization::CommodityLinearization, form::Formulation, M, N, rtrans, vtrans, ktrans, nv, nv_, na_, c, γc, γA, γt, γa1; sdtol=1e-10)
    # Linearization
    sumtx = custom_linearize_commodity_primal(model, linearization, primal(form), M, N, rtrans, vtrans, ktrans, nv, nv_, na_, c, γc, γA, γt, γa1)
    # Strong duality
    @constraint(model, sumtx <= unnormalized_objective_term(form) + sdtol)
    return sumtx
end

function custom_linearize_commodity_primal(model::Model, linearization::CommodityLinearization, primal::PrimalRepresentation, M, N, rtrans, vtrans, ktrans, nv, nv_, na_, c, γc, γA, γt, γa1)

    prob = problem(primal)
    parentprob = parent(prob)

    a1 = tolled_arcs(prob)
    k = index(prob)
    Amap = arcmap(prob)
	Vmap = used_nodes(prob)
	println("Vmap")
	
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

	
    b = NetPricing.sourcesink_vector(prob)			# Source sink vector 
    
 	bfull = expand_b(Vmap, nv, b)		# Source sink vector in full dimension
 	println("bfull")
 	γbfull_min, γbfull_avg, γbfull_max   = projection(vtrans, bfull)	# Projected Source sink vector full dimension
 	
 	# To verify if the projected problem is feasible
 	if γbfull_min == γbfull_max
 		γbfull = γbfull_min
 	else
 		γbfull = false
 	end
 	println(γbfull)
 	println("projection")
 	
 	if haskey(ktrans, γbfull)              # it is possible that the problem become infeasible in the transformed space
	
		k_ = ktrans[γbfull]               # retrieve the associated transformed problem
		println("k_")
		λ_ = rtrans.λvals[k_]                               # Reduce dual value in transformed space with the corresponding problem
		println("λ_")
		λ_full = expand_b(rtrans.Vmap[k_], nv_, λ_)         # Full dimension dual value in transformed space
		println("λ_full")
		γ_inv_λ_full = retroprojection(vtrans, value.(λ_full)) # Retroprojection of the dual value into the original space
		println("γ_inv_λ_full")
		x_ = rtrans.xvals[k_]
		println("x_")
		x_full = expand_b(rtrans.Amap[k_], na_, x_)       # optimal solution path in transformed problem 
		println("x_full")
		

		
		println("typeof, size λ_full      \t", typeof(λ_full), "\t "size(λ_full))
		println("typeof, size γA          \t", typeof(γA), "\t", size(γA))
		println("typeof, size c           \t", typeof(c), "\t", size(c))
		println("typeof, size γc          \t", typeof(γc), "\t", size(γc))
		println("typeof, size γt          \t", typeof(γt), "\t", size(γt))
		println("typeof, size x           \t", typeof(x), "\t", size(x))
		println("typeof, size b           \t", typeof(b), "\t", size(b))
		println("typeof, size γ_inv_λ_full\t", typeof(γ_inv_λ_full), "\t", size(γ_inv_λ_full))
		println("typeof, size γbfull      \t", typeof(γbfull), "\t", size(γbfull))

		prod1 = value.(γA' * λ_full)
		
		for i in γa1
			@constraint(model, prod1[i] ≤ γc[i] + γt[i])
			#@constraint(model, (γA' * λ_full)[i] ≤ γc[i]) # Is true for i not in γa1
		end
		println("γ(A)' * λ~ <= γ(c) + γ(t)")
		
		@constraint(model, c' * x + sumtx ≤ b' * γ_inv_λ_full)
		println("(c + t)' * x <= b' * γ^-1(λ~)")
		
		@constraint(model, (γc + γt)' * x_full ≤ γbfull' * λ_full)
		println("(γ(c) + γ(t))' * x~ <= γ(b)' * λ~")
	end
    

    # Run linearization-specific code
    linearize_commodity_extra(linearization, primal)

    return sumtx
end

