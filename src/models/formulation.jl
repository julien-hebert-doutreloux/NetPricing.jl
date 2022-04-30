#=
Abstract type of all formulations
A Formulation must define:
- problem(form): the overall problem (not the commodity problem)
- primal(form): primal representation
- dual(form): dual representation
- append!(model, form): add formulation to the model
- calculate_bigM(form): return a matrix whose row a corresponds to bigM candidates of arc a (order in a1)
- objective_term(form): term to add to the objective function
=#
abstract type Formulation end

primalobj(form::Formulation) = primalobj(primal(form))
dualobj(form::Formulation) = dualobj(dual(form))

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
