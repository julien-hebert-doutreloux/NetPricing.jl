## Structs
abstract type AbstractCommodityProblem <: AbstractProblem end

struct EmptyProblem <: AbstractCommodityProblem
    problem::Problem
    k::Int
end
    
struct UnprocessedProblem <: AbstractCommodityProblem
    problem::Problem
    k::Int
end

## Interfaces
function Base.parent(::AbstractCommodityProblem) end
function index(::AbstractCommodityProblem) end

orig(prob::AbstractCommodityProblem) = parent(prob).K[index(prob)].orig
dest(prob::AbstractCommodityProblem) = parent(prob).K[index(prob)].dest
demand(prob::AbstractCommodityProblem) = parent(prob).K[index(prob)].demand

arcmap(prob::AbstractCommodityProblem) = 1:length(arcs(parent(prob)))

paths(::AbstractCommodityProblem) = nothing

## Interface implementations
nodes(::EmptyProblem) = 0
arcs(::EmptyProblem) = ProblemArc[]
Base.parent(prob::EmptyProblem) = prob.problem
index(prob::EmptyProblem) = prob.k
used_nodes(::EmptyProblem) = Int[]
used_arcs(::EmptyProblem) = Int[]
paths(::EmptyProblem) = Vector{Int}[]

nodes(prob::UnprocessedProblem) = prob.problem.V
arcs(prob::UnprocessedProblem) = prob.problem.A
Base.parent(prob::UnprocessedProblem) = prob.problem
index(prob::UnprocessedProblem) = prob.k
used_nodes(prob::UnprocessedProblem) = collect(1:nodes(prob))
used_arcs(prob::UnprocessedProblem) = collect(1:length(arcs(prob)))

## Pretty print
function Base.show(io::IO, prob::EmptyProblem)
    print(io, "Empty problem for k = $(prob.k)")
end

function Base.show(io::IO, prob::UnprocessedProblem)
    print(io, "Unprocessed problem for k = $(prob.k)")
end