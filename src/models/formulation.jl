# Abstract type of all formulations
abstract type Formulation end

# Formulate: actualize all assigned formulations into model
function formulate(forms::Vector{<:Formulation}; silent=false, threads=nothing, kwargs...)
    model = Model(() -> Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), silent)
    set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)

    # Big M
    prob = problem(first(forms))
    M, N = calculate_bigM(prob)

    # Toll variables
    a1 = tolled_arcs(prob)
    a1dict = Dict(a => i for (i, a) in enumerate(a1))
    @variable(model, 0 ≤ t[a=a1], upper_bound = N[a1dict[a]])

    # Add formulations
    obj_terms = [append!(model, form, M, N; kwargs...) for form in forms]

    # Set objective function
    @objective(model, Max, sum(obj_terms))

    return model
end