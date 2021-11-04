module NetPricingTollBranch

using Printf
using JSON, Unmarshal
using Graphs, SimpleWeightedGraphs
using SparseArrays
using Hungarian
using JuMP, Gurobi
using DataStructures
# using GLMakie, GeometryBasics, Images
using ProgressMeter, Distributed
using BenchmarkTools, Random

import Base.Iterators: flatten, product
# import GLMakie.Axis

include("problems/abstractproblem.jl")
include("problems/problem.jl")
include("problems/preprocessedproblem.jl")
include("problems/tollsdict.jl")

include("misc/lazyenv.jl")

include("components/graph.jl")
include("components/graphtolls.jl")
include("components/bigm.jl")
include("components/model-components.jl")
include("components/graphmodel.jl")
include("components/fixedarcsolver.jl")
include("components/enumeration.jl")
include("components/preprocessing-path.jl")
include("components/preprocessing-spgm.jl")

# include("legacy/exactmodel.jl")
# include("legacy/upperbound.jl")
# include("legacy/lowerbound.jl")
# include("legacy/session.jl")
# include("legacy/logging.jl")
# include("legacy/branchandbound.jl")

# include("plots/plot.jl")
# include("plots/plot-multi.jl")

include("benchmark.jl")

export AbstractProblem, Problem, PreprocessedProblem, ProblemArc, Commodity
export read_problem, nodes, arcs, tolled_arcs, tollfree_arcs, srcdst_to_index, srcdst_to_cost, tolled_srcdst_to_index
export build_graph, shortest_path
export enumerate_bilevel_feasible
export preprocess_path, preprocess_spgm

# export BBSession, step!, solve!, solve
# export plot2d

export run_benchmark_fixedarcsolvers, run_benchmark_shortest_path, run_benchmark_enumeration

end # module
