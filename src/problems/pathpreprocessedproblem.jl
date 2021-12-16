## Path preprocessed problem
struct PathPreprocessedProblem <: AbstractPreprocessedProblem
    pprob::PreprocessedProblem
    paths::Vector{Vector{Int}}
    original_paths::Vector{Vector{Int}}
end

## Interfaces
nodes(preprob::PathPreprocessedProblem) = nodes(preprob.pprob)
arcs(preprob::PathPreprocessedProblem) = arcs(preprob.pprob)

Base.parent(prob::PathPreprocessedProblem) = parent(prob.pprob)
index(prob::PathPreprocessedProblem) = index(prob.pprob)
orig(prob::PathPreprocessedProblem) = orig(prob.pprob)
dest(prob::PathPreprocessedProblem) = dest(prob.pprob)

used_nodes(prob::PathPreprocessedProblem) = used_nodes(prob.pprob)
used_arcs(prob::PathPreprocessedProblem) = used_arcs(prob.pprob)

arcmap(prob::PathPreprocessedProblem) = arcmap(prob.pprob)

## Pretty print
function Base.show(io::IO, preprob::PathPreprocessedProblem)
    a1 = length(tolled_arcs(preprob))
    p = length(preprob.paths)
    print(io, "Path-processed problem for k = $(index(preprob)) with {$(nodes(preprob)) nodes, $(length(arcs(preprob))) arcs ($a1 tolled), $p paths}")
end