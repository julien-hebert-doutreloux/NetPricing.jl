## Structs
struct Commodity
    orig::Int
    dest::Int
    demand::Float64
end

struct Problem <: AbstractProblem
    V::Int
    A::Vector{ProblemArc}
    K::Vector{Commodity}
end

## Interfaces
nodes(prob::Problem) = prob.V
arcs(prob::Problem) = prob.A

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

## Sub-problem by commodity
subproblem_by_commodity(prob::Problem, k::Int) = Problem(prob.V, prob.A, [prob.K[k]])

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