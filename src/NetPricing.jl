module NetPricing

using Printf, Crayons, Crayons.Box, UnPack
using JSON, Unmarshal
using Graphs, SimpleWeightedGraphs
using SparseArrays
using Random, Distributions, StatsBase
using JuMP, Gurobi
using DataStructures, IterTools
using ProgressMeter, Distributed
using BenchmarkTools, Random

import Base.Iterators: flatten, product
import JuMP.Containers: DenseAxisArray

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
include("components/enumeration.jl")
include("components/preprocessing-path.jl")
include("components/preprocessing-spgm.jl")
include("components/preprocessing-light.jl")
include("components/preprocessing.jl")
include("components/conjugate-model.jl")

include("models/formulation.jl")
include("models/model-components.jl")
include("models/primal-representation.jl")
include("models/dual-representation.jl")

include("bigm/bigm-base.jl")
include("bigm/bigm-path.jl")

include("models/general-formulation.jl")
include("models/named-models.jl")
include("models/formulation-assignment.jl")

include("probgen/probgen.jl")
include("probgen/topologies.jl")

include("benchmark.jl")

export AbstractProblem, AbstractCommodityProblem, Problem, AbstractPreprocessedProblem, PreprocessedProblem, PathPreprocessedProblem, UnprocessedProblem, ProblemArc, Commodity
export read_problem, nodes, arcs, arcmap, paths, tolled_arcs, tollfree_arcs, srcdst_to_index, srcdst_to_cost, tolled_srcdst_to_index
export build_graph, shortest_path, get_path_cost, path_arcs, path_tolled_arcs
export enumerate_bilevel_feasible
export preprocess, preprocess_path, preprocess_spgm
export ConjugateModel, tvals, set_odpairs, set_demands, set_paths
export is_bilevel_feasible

export PrimalRepresentation, DualRepresentation, PrimalArc, PrimalPath, DualArc, DualPath
export Formulation, GeneralFormulation
export formulate!, assign, assign_breakpoint
export general_model, std_model, pastd_model, vf_model, pvf_model

export StandardFormulation, PathArcStandardFormulation, ValueFunctionFormulation, PathValueFunctionFormulation
export STDFormulation, PASTDFormulation, VFFormulation, PVFFormulation
export standard_model, path_arc_standard_model, value_function_model, path_value_function_model

export generate_problem, grid_graph

export consecutive_pairs

end # module
