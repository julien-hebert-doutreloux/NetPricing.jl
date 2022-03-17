module NetPricing

using Printf, Crayons, Crayons.Box
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
include("problems/abstractcommodityproblem.jl")
include("problems/preprocessedproblem.jl")
include("problems/pathpreprocessedproblem.jl")
include("problems/tollsdict.jl")

include("misc/lazyenv.jl")
include("misc/utils.jl")

include("components/graph.jl")
include("components/graphtolls.jl")
include("components/graphdists.jl")
include("components/graphmodel.jl")
include("components/fixedarcsolver.jl")
include("components/bigm.jl")
include("components/smallm.jl")
include("components/enumeration.jl")
include("components/preprocessing-path.jl")
include("components/preprocessing-spgm.jl")
include("components/preprocessing-light.jl")
include("components/preprocessing.jl")

include("models/model-components.jl")
include("models/refs-macros.jl")
include("models/standard-model.jl")

# include("plots/plot.jl")
# include("plots/plot-multi.jl")

include("benchmark.jl")

export AbstractProblem, AbstractCommodityProblem, Problem, AbstractPreprocessedProblem, PreprocessedProblem, PathPreprocessedProblem, UnprocessedProblem, ProblemArc, Commodity
export read_problem, nodes, arcs, arcmap, paths, tolled_arcs, tollfree_arcs, srcdst_to_index, srcdst_to_cost, tolled_srcdst_to_index
export build_graph, shortest_path, get_path_cost
export enumerate_bilevel_feasible
export preprocess, preprocess_path, preprocess_spgm
export standard_model

# export plot2d

export run_benchmark_fixedarcsolvers, run_benchmark_shortest_path, run_benchmark_enumeration

end # module
