#=
Abstract type of all formulations
A Formulation must define:
- problem(form): the overall problem (not the commodity problem)
- primal(form): primal representation
- dual(form): dual representation
- append!(model, form): add formulation to the model
- calculate_bigM(form): return a matrix whose row a corresponds to bigM candidates of arc a (order in a1)
- objective_term(form): term to add to the objective function
- unnormalized_objective_term(form): objective term before multiply by demand
=#
abstract type Formulation end

abstract type AbstractLinearization end

# Formulate: actualize all assigned formulations into model
function formulate!(forms::Vector{<:Formulation}, linearization::AbstractLinearization; silent=false, threads=nothing, sdtol=1e-10, kwargs...)
    model = Model(() -> Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), silent)
    set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)

    prob = problem(first(forms))

    # Big M
    Ms = [calculate_bigM(form, threads=threads) for form in forms]
    N = max.(collect.(maximum.(Ms, dims=2))...)

    # Toll variables
    a1 = tolled_arcs(prob)
    a1dict = Dict(a => i for (i, a) in enumerate(a1))
    @variable(model, 0 ≤ t[a=a1], upper_bound = N[a1dict[a]])

    # Add formulations
    for form in forms
        append!(model, form; kwargs...)
    end

    # Linearize (and strong duality)
    linearize!(model, linearization, forms, Ms, N; sdtol=sdtol)

    # Set objective function
    @objective(model, Max, sum(objective_term.(forms)))

    return model, forms
end



# Formulate: actualize all assigned formulations into model (CUSTOM)
# function custom_formulate!(forms::Vector{<:Formulation}, linearization::AbstractLinearization, Ms, NL, option; silent=false, threads=nothing, sdtol=1e-10, kwargs...)
#     model = Model(() -> Gurobi.Optimizer(current_env())) # NetPricing.
#     set_optimizer_attribute(model, MOI.Silent(), silent)
#     set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)

#     prob = problem(first(forms))

#     # Big M
#     #Ms = [calculate_bigM(form, threads=threads) for form in forms]
#     NU = max.(collect.(maximum.(Ms, dims=2))...)

#     # Toll variables
#     a1 = tolled_arcs(prob)
#     a1dict = Dict(a => i for (i, a) in enumerate(a1))
    
#     if option == 0
# 		# normal
# 		@variable(model, 0 ≤ t[a=a1], upper_bound = NU[a1dict[a]])
# 		N = NU
# 		#println("option 0")
#     elseif option == 1
#     	# shortest path
#     	@variable(model, t[a=a1], lower_bound = NL[a1dict[a]], upper_bound = NL[a1dict[a]])
# 		N = NU
# 		#println("option 1")
# 	elseif option == 2
# 		# lower bound
# 		@variable(model, t[a=a1], lower_bound = NL[a1dict[a]], upper_bound = NU[a1dict[a]])
# 		N = NU
# 		#println("option 2")
# 	elseif option == 3 
# 		# upper bound
# 		@variable(model, 0 ≤ t[a=a1], upper_bound = NL[a1dict[a]])
# 		N = NL
# 		#println("option 3")
# 	elseif option == 4
# 		# comprehensive lower bound
# 		L = Int.(NL .< NU) .* NL
# 		@variable(model, t[a=a1], lower_bound = L[a1dict[a]], upper_bound = NU[a1dict[a]])
# 		N = NU
# 		#println("option 4")
# 	elseif option == 5
# 		# comprehensive upper bound
# 		U = Int.(NL .< NU) .* NL + Int.(NU .<= NL) .* NU
# 		@variable(model, 0 ≤ t[a=a1], upper_bound = U[a1dict[a]])
# 		N = U
# 		#println("option 5")
# 	end
	
#     # Add formulations
#     for form in forms
#         append!(model, form; kwargs...)
#     end

#     # Linearize (and strong duality)
#     linearize!(model, linearization, forms, Ms, N; sdtol=sdtol) #NetPricing.

#     # Set objective function
#     @objective(model, Max, sum(objective_term.(forms))) # NetPricing.

#     return model, forms
# end


function custom_formulate!(forms::Vector{<:Formulation}, linearization::AbstractLinearization, Ms, NL, option; rtrans=nothing, trans=nothing, silent=false, threads=nothing, sdtol=1e-10, kwargs...)
    model = Model(() -> Gurobi.Optimizer(current_env())) # NetPricing.
    set_optimizer_attribute(model, MOI.Silent(), silent)
    set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)
    
    prob = problem(first(forms))

    # Big M
    #Ms = [calculate_bigM(form, threads=threads) for form in forms]
    NU = max.(collect.(maximum.(Ms, dims=2))...)

    # Toll variables
    a1 = tolled_arcs(prob)
    a1dict = Dict(a => i for (i, a) in enumerate(a1))
    
    if option == 0 || option == 6
		# normal
		@variable(model, 0 ≤ t[a=a1], upper_bound = NU[a1dict[a]])
		N = NU
		#println("option 0")
    elseif option == 1
    	# shortest path
    	@variable(model, t[a=a1], lower_bound = NL[a1dict[a]], upper_bound = NL[a1dict[a]])
		N = NU
		#println("option 1")
	elseif option == 2
		# lower bound
		@variable(model, t[a=a1], lower_bound = NL[a1dict[a]], upper_bound = NU[a1dict[a]])
		N = NU
		#println("option 2")
	elseif option == 3 
		# upper bound
		@variable(model, 0 ≤ t[a=a1], upper_bound = NL[a1dict[a]])
		N = NL
		#println("option 3")
	elseif option == 4
		# comprehensive lower bound
		L = Int.(NL .< NU) .* NL
		@variable(model, t[a=a1], lower_bound = L[a1dict[a]], upper_bound = NU[a1dict[a]])
		N = NU
		#println("option 4")
	elseif option == 5
		# comprehensive upper bound
		U = Int.(NL .< NU) .* NL + Int.(NU .<= NL) .* NU
		@variable(model, 0 ≤ t[a=a1], upper_bound = U[a1dict[a]])
		N = U
		#println("option 5")
	end
	if (option == 6) && (rtrans != nothing) && (trans != nothing)
	
		vtrans = trans["V"]# Vertex transformation (before:after)
		etrans = trans["A"]# Edge transformation (before:after)
		
		c  = cost_vector(prob)      # constant cost vector from the original problem
		γc = projection(etrans, c)  # Projection of c in transformed space
		println("γc")
		
		
		γA = trans["M_"]    # incidence matrix in transformed space
		γA = hcat(γA...)'   # convert in a matrix (vertex, edge)

		nv_ = size(γA)[1]   # number of node in the transformed space
		na_ = size(γA)[2]   # number of arcs in the transformed space
		###############
		

		ktrans = Dict()     # mapping between full b~ and its associated index k (b~:after)
		# Only works for result in transformed space (when x is in the id)
		for (k, v) in rtrans.b
			VV = rtrans.Vmap[k]
			bfull = expand_b(VV, nv_, v) 
			ktrans[bfull] = k
		end
		println("ktrans")

		γa1 = []         # tolled edge index in the transformed space
		γa1dict = Dict() # mapping entre indices et classes d'équivalences d'indice
		for e in trans["RA"]
			if !isempty(intersect(e, a1))
				println(e[1])
				i = etrans[e[1]]
				append!(γa1, i)
				γa1dict[i] = e
			end
		end
		println("γa1dict")
		
		# variable artificiel γt
		@variable(model, γt[a=γa1], base_name="γt") # pas besoin de borner voir les contraintes
		for (k, v) in γa1dict
			@constraint(model, γt[k] == mean(t[v]))
		end
				
		# γ(A), λ~, γ(c), c, b, γ^-1(λ~), x~, γ(b) sont des constantes
		# γ(A)' * λ~ <= γ(c) + γ(t) need linearize_commodity_primal_custom
		# (c + t)' * x <= b' * γ^-1(λ~) need linearize_commodity_primal_custom
		# (γ(c) + γ(t))' * x~ <= γ(b)' * λ~ need linearize_commodity_primal_custom
	end
	
    # Add formulations
    for form in forms
        append!(model, form; kwargs...)
    end

	if option==6
		custom_linearize!(model, linearization, forms, Ms, N, rtrans, vtrans, ktrans, nv_,na_, c, γc, γA, γt; sdtol=sdtol)
	else
		# Linearize (and strong duality)
		linearize!(model, linearization, forms, Ms, N; sdtol=sdtol) #NetPricing.
	end
	
	
    # Set objective function
    @objective(model, Max, sum(objective_term.(forms))) # NetPricing.

    return model, forms
end

