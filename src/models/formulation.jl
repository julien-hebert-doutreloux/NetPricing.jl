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
function custom_formulate!(forms::Vector{<:Formulation}, linearization::AbstractLinearization, Ms, NL, option; silent=false, threads=nothing, sdtol=1e-10, kwargs...)
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
    
    if option == 0
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
	
    # Add formulations
    for form in forms
        append!(model, form; kwargs...)
    end

    # Linearize (and strong duality)
    linearize!(model, linearization, forms, Ms, N; sdtol=sdtol) #NetPricing.

    # Set objective function
    @objective(model, Max, sum(objective_term.(forms))) # NetPricing.

    return model, forms
end
