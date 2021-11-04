## Structs
struct PreprocessedProblem <: AbstractProblem
    problem::Problem
    V::Int
    A::Vector{ProblemArc}
    Vmap::Vector{Int}
    Amap::Vector{Int}
    Vrevmap::Vector{Int}
    Arevmap::Vector{Int}
    k::Int
    orig::Int
    dest::Int
end

## Interfaces
nodes(preprob::PreprocessedProblem) = preprob.V
arcs(preprob::PreprocessedProblem) = preprob.A

## Reversed map
function revmap(map::AbstractVector{Int}, maxlength::Int)
    rmap = zeros(Int, maxlength)
    for (i, v) in enumerate(map)
        rmap[v] = i
    end
    return rmap
end

## Pretty print
function Base.show(io::IO, preprob::PreprocessedProblem)
    a1 = length(tolled_arcs(preprob))
    print(io, "PreprocessedProblcomem for k = $(preprob.k) with {$(preprob.V) nodes, $(length(preprob.A)) arcs ($a1 tolled)}")
end