# NetPricing <!-- omit from toc -->

Julia library for the Network Pricing Problem (NPP). For more details, read our papers
\[[1](#readme-ref1), [2](#readme-ref2)\].

## Tables of Contents <!-- omit from toc -->

- [Installation](#installation)
- [Basic Usage](#basic-usage)
- [Example of a Problem File](#example-of-a-problem-file)
- [Other Features](#other-features)
  - [Shortest Path Solver](#shortest-path-solver)
  - [Forward and Conjugate Models](#forward-and-conjugate-models)
  - [Problem Generation](#problem-generation)
  - [Bilevel Feasibility Test](#bilevel-feasibility-test)
  - [Strong Bilevel Feasibility Test](#strong-bilevel-feasibility-test)
- [Problem Instance Sets](#problem-instance-sets)
- [References](#references)

## Installation

The package is not available on the registry. Please clone the repository and
add it as a [local package](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-a-local-package). You also need to add [JuMP](https://jump.dev/JuMP.jl/stable/) and [Gurobi](https://github.com/jump-dev/Gurobi.jl).

## Basic Usage

```julia
using NetPricing, JuMP, Gurobi

# Import a problem from a file
file = "./problems/comm-heavy/i120-01.json"
prob = read_problem(file)

# Preprocess the problem for each commodity
pprobs = preprocess(prob, maxpaths = 1000)

# Create a model
# model, forms = std_model(pprobs)
model, forms = pastd_model(pprobs)
# model, forms = vf_model(pprobs)
# model, forms = pvf_model(pprobs)

# Add strong bilevel feasibility (optional, only available for pastd and pvf models)
add_strong_bf_cuts(model, forms, maxpairs=10000, commpairs=100)

# Solve the model
optimize!(model)

# Extract the result
tvals = value.(model[:t])                   # The prices t

k = 1
primal_repr = primal(forms[k])              # Primal representation
dual_repr = NetPricing.dual(forms[k])       # Dual representation
prob_k = problem(primal_repr)               # Preprocessed problem of forms[k]

primal_obj = value(primal_repr.primalobj)   # Primal objective: c' x[k]
dual_obj = value(dual_repr.dualobj)         # Dual objective: b' λ[k]

xvals = value.(primal_repr.x)               # Arc selections x[k]
zvals = value.(primal_repr.z)               # Path selections z[k] (only for primal-path)
λvals = value.(dual_repr.λ)                 # Dual prices λ[k] (only for dual-arc)

# Maps from indices of the preprocessed problem to those of the base problem
Amap = used_arcs(prob_k)
Vmap = used_nodes(prob_k)

# Primal and dual objectives must satisfy: primal_obj + t' x[k] == dual_obj
a1_k = tolled_arcs(prob_k)                  # List of tolled arcs of prob_k
@assert(dual_obj - primal_obj ≈ sum(tvals[Amap[a]] .* xvals[a] for a in a1_k))
```

**Note**: `x` and `λ` are indexed based on the preprocessed problem. To map back to the base problem, use `Amap` and `Vmap`. In particular, `x[a]` corresponds to arc `Amap[a]` in the base problem, and `λ[v]` corresponds `Vmap[v]`. If the primal representation is primal-path, then the path corresponding to `z[p]` is `paths(prob_k)[p]` which is a list of nodes indexed based on the base problem.

## Example of a Problem File
This is a problem file (JSON format) describing an NPP instance with 3 nodes (field `V`),
3 arcs (field `A`) consisting of 1 tolled arc (arc `1-3`) and 2 toll-free arcs (arcs `1-2` and `2-3`), and 1 commodity (field `K`) from node `1` to node `3`. Note that Julia uses 1-based numbering, ie. array indices start from 1.
```json
{
    "problem": {
        "V": 3,
        "A": [{
            "src": 1,
            "dst": 3,
            "cost": 1.0,
            "toll": true
        }, {
            "src": 1,
            "dst": 2,
            "cost": 2.0,
            "toll": false
        }, {
            "src": 2,
            "dst": 3,
            "cost": 3.0,
            "toll": false
        }],
        "K": [{
            "orig": 1,
            "dest": 3,
            "demand": 1.0
        }]
    }
}
```

## Other Features
### Shortest Path Solver
```julia
file = "./problems/comm-heavy/i120-01.json"
prob = read_problem(file)

# List of indices of tolled arcs
a1 = tolled_arcs(prob)

# Shortest path
graph = NetPricing.build_graph(prob)
orig, dest = 1, 36
p_nodes, p_cost = shortest_path(graph, orig, dest)

# List of arcs (tolled_arcs) along the path p_nodes
p_arcs = collect(path_arcs(p_nodes, prob))
p_tolled_arcs = path_tolled_arcs(p_nodes, prob)

# Shortest path with tolls
NetPricing.set_tolls!(graph, prob, Dict(106 => 10., 109 => 20.))
p_nodes, p_cost = shortest_path(graph, orig, dest)

# Toll-free path (t = infty)
NetPricing.disable_tolls!(graph, prob)
p_nodes, p_cost = shortest_path(graph, orig, dest)

# Null-toll path (t = 0)
NetPricing.reset_tolls!(graph, prob)
p_nodes, p_cost = shortest_path(graph, orig, dest)
```

### Forward and Conjugate Models
```julia
file = "./problems/comm-heavy/i120-01.json"
prob = read_problem(file)
a1 = tolled_arcs(prob)
pprobs = preprocess(prob)

# Forward model (solve for w given t)
fsolver = ForwardHybridSolver(pprobs)
t_dict = Dict(a1 .=> 1.)
ft, w_vec = solve(fsolver, t_dict)      # Value of f(t) and w in Vector form
w_dict = Dict(a1 .=> w_vec)             # w in Dict form

# Conjugate model (solve for t given w)
csolver = ConjugateDynamicHybridModel(pprobs)
w_dict = Dict(a1 .=> 1.)
gw, t_vec = solve(csolver, w_dict)      # Value of g(w) and t in Vector form
t_dict = Dict(a1 .=> t_vec)             # t in Dict form
```

### Problem Generation
```julia
# Choose a graph structure
graph = grid_graph(6, 6)

# Generate a random problem
num_commodities = 120
prob = generate_problem(graph, num_commodities)

# Create a problem from scratch
V = 3
A = [
    ProblemArc(1, 3, 1., true),     # src, dst, cost, toll
    ProblemArc(1, 2, 2., false),
    ProblemArc(2, 3, 3., false)]
K = [
    Commodity(1, 3, 1.)             # orig, dest, demand
    ]
prob = Problem(V, A, K)
```
### Bilevel Feasibility Test
```julia
file = "./problems/partialparallel3.json"
prob = read_problem(file)

# Bilevel feasibility test (for 2 commodities)
nk = 2
tester = ConjugatePrimalTester(prob, nk)
p, q1, q2 = [1,11,12,6], [3,11,12,8], [3,13,12,8]

is_bilevel_feasible(tester, [p, q1])      # true
is_bilevel_feasible(tester, [p, q2])      # false
```

### Strong Bilevel Feasibility Test
This is a demonstration of Example 9 in \[[2](#readme-ref2)\].
```julia
file = "./test/problems/bftest-triplecomm.json"
prob = read_problem(file)

tester2 = ConjugatePrimalTester(prob, 2)
tester3 = ConjugatePrimalTester(prob, 3)
p, q, r = [1,7,8,4], [2,5], [3,9,10,6]

is_strongly_bilevel_feasible(tester2, [p, q])       # true
is_strongly_bilevel_feasible(tester2, [p, r])       # true
is_strongly_bilevel_feasible(tester2, [q, r])       # true

is_strongly_bilevel_feasible(tester3, [p, q, r])    # false
```

## Problem Instance Sets
The problem instances tested in \[[1](#readme-ref1)\] are in the folder `./problems/paper/`. There are 4 classes: class `g` (`5x12` grid), class `h` (`12x12` grid), class `d` (`144`-node Delaunay), and class `v` (`144`-node Voronoi).

Instances tested in \[[2](#readme-ref2)\] are in the folder `./problems/progressive-2/`. The grid size increases gradually from class `i` (`6x6`) to class `o` (`12x12`).

## References
<a id="readme-ref1"></a> \[1\] Quang Minh Bui, Bernard Gendron, Margarida Carvalho. A Catalog of Formulations for the Network Pricing Problem. INFORMS Journal on Computing, 34(5):2658–2674, 2022. ([arXiv](https://arxiv.org/abs/2106.03887))

<a id="readme-ref2"></a> \[2\] Quang Minh Bui, Margarida Carvalho, José Neto. Asymmetry in the Complexity of the Multi-Commodity Network Pricing Problem. ([arXiv](https://arxiv.org/abs/2212.10673))