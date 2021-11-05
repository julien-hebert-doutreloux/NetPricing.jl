## Structs
abstract type AbstractProblem end

struct ProblemArc
    src::Int
    dst::Int
    cost::Float64
    toll::Bool
end

## Interfaces
function nodes(::AbstractProblem) end
function arcs(::AbstractProblem) end

## Queries
tolled_arcs(prob::AbstractProblem) = findall(a -> a.toll, arcs(prob))
tollfree_arcs(prob::AbstractProblem) = findall(a -> !a.toll, arcs(prob))

srcdst_to_index(prob::AbstractProblem) = Dict((a.src, a.dst) => i for (i, a) in enumerate(arcs(prob)))
srcdst_to_cost(prob::AbstractProblem) = Dict((a.src, a.dst) => a.cost for a in arcs(prob))

tolled_srcdst_to_index(prob::AbstractProblem) = Dict((a.src, a.dst) => i for (i, a) in enumerate(arcs(prob)) if a.toll)
