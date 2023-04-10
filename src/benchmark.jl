using BenchmarkTools
using Random

function run_benchmark_shortest_path(probfile="problems/paper/d30-01.json")
    prob = read_problem(probfile)
    graph = NetPricing.build_graph(prob)
    
    println("Built-in Dijkstra")
    tNormal = @benchmark NetPricing.shortest_path_old($graph, c.orig, c.dest) setup=(c = rand($prob.K))
    display(tNormal)
    println()

    println("Minimal Dijkstra")
    tMinimal = @benchmark NetPricing.shortest_path($graph, c.orig, c.dest) setup=(c = rand($prob.K))
    display(tMinimal)
    println()

    println(judge(median(tMinimal), median(tNormal)))
end

function run_benchmark_enumeration(probfile="problems/paper/d30-01.json", numpaths=10000; kwargs...)
    prob = read_problem(probfile)
    graph = NetPricing.build_graph(prob)
    totalnumpaths = 0

    time = @elapsed begin
        for (k, comm) in enumerate(prob.K)
            paths = enumerate_bilevel_feasible(graph, comm.orig, comm.dest, prob, numpaths; kwargs...)
            totalnumpaths += length(paths)

            @printf "%6d " length(paths)
            (k % 10 == 0) && println()
        end
    end

    println()
    @printf "Total: %-d paths\n" totalnumpaths
    @printf "Time:  %-.2f s\n" time
    @printf "Speed: %-.0f paths/s\n" (totalnumpaths / time)
end

function run_benchmark_conjugate_solver(probfile="problems/paper/d30-01.json"; threads=nothing)
    prob = read_problem(probfile)
    a1 = tolled_arcs(prob)
    a1set = BitSet(a1)

    println("Preprocessing...")
    pprobs = preprocess(prob, maxpaths=10000)
    max_w = Dict(a1 .=> 0.)
    for pprob in pprobs
        Amap = used_arcs(pprob)
        for a in Amap ∩ a1set
            max_w[a] += demand(pprob)
        end
    end
    println()

    function make_demands()
        return Dict(a => rand() * max_w[a] for a in a1)
    end
    
    println("Linear Solver")
    sLinear = NetPricing.ConjugateLinearModel(prob; threads)
    tLinear = @benchmark solve($sLinear, d) setup=(d = $make_demands())
    display(tLinear)
    println()

    println("Dynamic Linear Solver")
    sDynamic = NetPricing.ConjugateDynamicLinearModel(prob; threads)
    tDynamic = @benchmark solve($sDynamic, d) setup=(d = $make_demands())
    display(tDynamic)
    println()

    for maxpaths in [0, 100, 1000, 10000]
        println("Dynamic Hybrid Solver - $maxpaths paths")
        sDynamicHybrid = NetPricing.ConjugateDynamicHybridModel(pprobs; maxpaths, threads)
        tDynamicHybrid = @benchmark solve($sDynamicHybrid, d) setup=(d = $make_demands())
        display(tDynamicHybrid)
        println()
    end
    
    println("Preprocessed DualArc Solver")
    sPreprocArc = NetPricing.ConjugatePreprocessedModel(DualArc, pprobs; threads)
    tPreprocArc = @benchmark solve($sPreprocArc, d) setup=(d = $make_demands())
    display(tPreprocArc)
    println()

    println("Preprocessed DualPath Solver")
    sPreprocPath = NetPricing.ConjugatePreprocessedModel(DualPath, pprobs; threads)
    tPreprocPath = @benchmark solve($sPreprocPath, d) setup=(d = $make_demands())
    display(tPreprocPath)
    println()

    # println("KKT Solver")
    # sKKT = NetPricing.ConjugateKKTModel(prob)
    # tKKT = @benchmark solve($sKKT, d) setup=(d = $make_demands())
    # display(tKKT)
    # println()
end

function run_benchmark_forward_solver(probfile="problems/paper/d30-01.json")
    prob = read_problem(probfile)
    a1 = tolled_arcs(prob)

    println("Preprocessing...")
    pprobs = preprocess(prob, maxpaths=10000)
    filter!(p -> !(p isa EmptyProblem), pprobs)

    println("Calculate big-M...")
    Ms = [prob isa PathPreprocessedProblem ? calculate_bigM_difference(prob) : calculate_bigM(prob) for prob in pprobs]
    N = Dict(a1 .=> max.(collect.(maximum.(Ms, dims=2))...))
    println()

    function make_tolls()
        return Dict(a => rand() * N[a] for a in a1)
    end

    for maxpaths in [0, 10, 100, 1000, 2000, 5000, 10000]
        println("Hybrid Solver - $maxpaths paths")
        solver = NetPricing.ForwardHybridSolver(pprobs; maxpaths)
        time = @benchmark solve($solver, t) setup=(t = $make_tolls())
        display(time)
        println()
    end
end
