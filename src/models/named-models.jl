std_model(probs; kwargs...) = general_model(probs, primal_arc, dual_arc; kwargs...)
const standard_model = std_model

pastd_model(probs; kwargs...) = general_model(probs, primal_path, dual_arc; kwargs...)
const patharc_standard_model = pastd_model