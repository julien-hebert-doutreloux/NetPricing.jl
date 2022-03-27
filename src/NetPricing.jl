module NetPricing

using Printf, Crayons, Crayons.Box, UnPack
using JSON, Unmarshal
using Graphs, SimpleWeightedGraphs
using SparseArrays
using Hungarian
using Random, Distributions, StatsBase
using JuMP, Gurobi
import JuMP.Containers: DenseAxisArray
using DataStructures, IterTools
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
include("misc/consecutive-pairs.jl")

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
include("components/inverse-model.jl")

include("bigm-path/bigm-path.jl")

include("models/model-components.jl")
include("models/refs-macros.jl")
include("models/formulation.jl")
include("models/general-formulation.jl")
include("models/formulation-assignment.jl")
include("models/primal-representation.jl")
include("models/dual-representation.jl")
include("models/named-models.jl")

include("probgen/probgen.jl")
include("probgen/topologies.jl")

# include("plots/plot.jl")
# include("plots/plot-multi.jl")

include("benchmark.jl")

export AbstractProblem, AbstractCommodityProblem, Problem, AbstractPreprocessedProblem, PreprocessedProblem, PathPreprocessedProblem, UnprocessedProblem, ProblemArc, Commodity
export read_problem, nodes, arcs, arcmap, paths, tolled_arcs, tollfree_arcs, srcdst_to_index, srcdst_to_cost, tolled_srcdst_to_index
export build_graph, shortest_path, get_path_cost, path_arcs, path_tolled_arcs
export enumerate_bilevel_feasible
export preprocess, preprocess_path, preprocess_spgm
export inverse_model, set_inverse_model_odpairs, set_inverse_model_paths

export PrimalRepresentation, DualRepresentation, PrimalArc, PrimalPath, DualArc, DualPath
export Formulation, GeneralFormulation
export StandardFormulation, PathArcStandardFormulation, ValueFunctionFormulation, PathValueFunctionFormulation
export STDFormulation, PASTDFormulation, VFFormulation, PVFFormulation
export formulate, assign, assign_breakpoint
export general_model, std_model, pastd_model, vf_model, pvf_model
export standard_model, path_arc_standard_model, value_function_model, path_value_function_model

export generate_problem, generate_cost
export grid_graph

export consecutive_pairs

# export plot2d

export run_benchmark_fixedarcsolvers, run_benchmark_shortest_path, run_benchmark_enumeration

end # module
