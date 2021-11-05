## Structs
struct PreprocessedProblem <: AbstractCommodityProblem
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

Base.parent(prob::PreprocessedProblem) = prob.problem
index(prob::PreprocessedProblem) = prob.k
orig(prob::PreprocessedProblem) = prob.orig
dest(prob::PreprocessedProblem) = prob.dest

arcmap(prob::PreprocessedProblem) = prob.Amap

## Reversed map
function revmap(map::AbstractVector{Int}, maxlength::Int)
    rmap = zeros(Int, maxlength)
    for (i, v) in enumerate(map)
        (v in 1:maxlength) && (rmap[v] = i)
    end
    return rmap
end

## Pretty print
function Base.show(io::IO, preprob::PreprocessedProblem)
    a1 = length(tolled_arcs(preprob))
    print(io, "Processed problem for k = $(preprob.k) with {$(preprob.V) nodes, $(length(preprob.A)) arcs ($a1 tolled)}")
end