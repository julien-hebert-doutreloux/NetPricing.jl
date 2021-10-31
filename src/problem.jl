## Structs
struct ProblemArc
    src::Int
    dst::Int
    cost::Float64
    toll::Bool
end

struct Commodity
    orig::Int
    dest::Int
    demand::Float64
end

struct Problem
    V::Int
    A::Vector{ProblemArc}
    K::Vector{Commodity}
end

## Write to file
function Base.write(filename::AbstractString, prob::Problem)
    open(filename, "w") do f
        JSON.print(f, Dict(:problem => prob))
    end
end

## Read from file
function read_problem(filename::AbstractString)
    local prob
    open(filename, "r") do f
        str = read(f, String)
        prob = unmarshal(Problem, JSON.parse(str)["problem"])
    end
    return prob
end

## Queries
tolled_arcs(prob::Problem) = [i for i in 1:length(prob.A) if prob.A[i].toll]
tollfree_arcs(prob::Problem) = [i for i in 1:length(prob.A) if !prob.A[i].toll]

## Sub-problem by commodity
subproblem_by_commodity(prob::Problem, k::Int) = Problem(prob.V, prob.A, [prob.K[k]])

## Dicts
srcdst_to_index(prob::Problem) = Dict((a.src, a.dst) => i for (i, a) in enumerate(prob.A))

## Pretty print
function Base.show(io::IO, prob::Problem)
    a1 = length(tolled_arcs(prob))
    print(io, "Problem with {$(prob.V) nodes, $(length(prob.A)) arcs ($a1 tolled), $(length(prob.K)) commodities}")
end

function Base.show(io::IO, arc::ProblemArc)
    tollstr = arc.toll ? "Tolled" : "Toll-free"
    print(io, "$tollstr arc $(arc.src) => $(arc.dst) with cost $(arc.cost)")
end

function Base.show(io::IO, comm::Commodity)
    print(io, "Commodity $(comm.orig) => $(comm.dest) with demand $(comm.demand)")
end