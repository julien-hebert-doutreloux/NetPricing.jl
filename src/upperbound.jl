## Build model
function make_upperbound_model(prob::Problem; kwargs...)
    model = make_exact_model(prob, kwargs...)
    unset_binary.(model[:x][:,model[:a1]])
    return model
end
