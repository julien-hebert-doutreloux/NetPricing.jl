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
using Parameters: @with_kw
using VoronoiDelaunay

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
include("misc/binary-search.jl")

include("components/graph.jl")
include("components/graphtolls.jl")
include("components/graphdists.jl")
include("components/graphmodel.jl")
include("components/enumeration.jl")
include("components/preprocessing-path.jl")
include("components/preprocessing-spgm.jl")
include("components/preprocessing-light.jl")
include("components/preprocessing.jl")

include("models/formulation.jl")
include("models/model-components.jl")
include("models/primal-representation.jl")
include("models/dual-representation.jl")

include("forward-solvers/abstract.jl")
include("forward-solvers/abstract-commodity.jl")
include("forward-solvers/graph-solver.jl")
include("forward-solvers/path-solver.jl")
include("forward-solvers/hybrid-solver.jl")

include("conjugate-solvers/abstract.jl")
include("conjugate-solvers/linear-model.jl")
include("conjugate-solvers/dynamic-linear-model.jl")
include("conjugate-solvers/dynamic-hybrid-model.jl")
include("conjugate-solvers/kkt-model.jl")
include("conjugate-solvers/preprocessed-model.jl")

include("bf-tests/abstract.jl")
include("bf-tests/conjugate-tester.jl")
include("bf-tests/conjugate-kkt-tester.jl")
include("bf-tests/conjugate-primal-tester.jl")

include("bigm/bigm-base.jl")
include("bigm/bigm-path.jl")
include("bigm/bigm-difference.jl")

include("models/general-formulation.jl")
include("models/linearizations/commodity-linearization.jl")
include("models/linearizations/binary-linearization.jl")
include("models/named-models.jl")
include("models/formulation-assignment.jl")

include("strong-bf-cuts/biclique-cover.jl")
include("strong-bf-cuts/cut-generation.jl")
include("strong-bf-cuts/model-augmentation.jl")

include("probgen/probgen-args.jl")
include("probgen/probgen.jl")
include("probgen/topologies.jl")
include("probgen/progressive.jl")

include("misc/default.jl")

include("benchmark.jl")

include("misc/gamma.jl")
export expand_b, retroprojection, projection



export AbstractProblem, AbstractCommodityProblem, Problem, AbstractPreprocessedProblem, PreprocessedProblem, PathPreprocessedProblem, UnprocessedProblem, ProblemArc, Commodity
export read_problem, nodes, arcs, arcmap, paths, tolled_arcs, tollfree_arcs, srcdst_to_index, srcdst_to_cost, tolled_srcdst_to_index
export used_nodes, used_arcs
export build_graph, shortest_path, get_path_cost, path_arcs, path_tolled_arcs
export enumerate_bilevel_feasible
export preprocess, preprocess_path, preprocess_spgm

export AbstractForwardSolver, set_tolls, wvals, solve
export ForwardHybridSolver

export AbstractConjugateSolver, set_odpairs, set_demands, set_paths, tvals, solve
export ConjugateLinearModel, ConjugateDynamicLinearModel, ConjugateDynamicHybridModel, ConjugateKKTModel, ConjugatePreprocessedModel

export AbstractBFTester, AbstractStrongBFTester, is_bilevel_feasible, is_strongly_bilevel_feasible
export ConjugateBFTester, ConjugateKKTTester, ConjugatePrimalTester

export PrimalRepresentation, DualRepresentation, PrimalArc, PrimalPath, DualArc, DualPath
export AbstractLinearization, CommodityLinearization, ArcLinearization, PathLinearization, BinaryLinearization
export Formulation, GeneralFormulation
export formulate!, assign, assign_breakpoint
export general_model, std_model, pastd_model, vf_model, pvf_model
export problem, primal, dual

export StandardFormulation, PathArcStandardFormulation, ValueFunctionFormulation, PathValueFunctionFormulation
export STDFormulation, PASTDFormulation, VFFormulation, PVFFormulation
export standard_model, path_arc_standard_model, value_function_model, path_value_function_model

export add_strong_bf_cuts

export generate_problem, grid_graph, delaunay_graph, voronoi_graph
export generate_progressive_grid, generate_progressive_delaunay

export consecutive_pairs

end # module
