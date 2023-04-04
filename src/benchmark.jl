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

function run_benchmark_conjugate_solver(probfile="problems/paper/d30-01.json")
    prob = read_problem(probfile)
    a1 = tolled_arcs(prob)
    a1set = BitSet(a1)

    println("Preprocessing...")
    pprobs = preprocess(prob, maxpaths=10000)
    maxk = Dict(a1 .=> 0)
    for pprob in pprobs
        Amap = arcmap(pprob)
        for a in Amap ∩ a1set
            maxk[a] += 1
        end
    end
    println()

    function make_demands()
        return Dict(a => rand(0:maxk[a]) for a in a1)
    end
    
    println("Linear Solver")
    sLinear = NetPricing.ConjugateLinearModel(prob)
    tLinear = @benchmark solve($sLinear, d) setup=(d = $make_demands())
    display(tLinear)
    println()

    println("Dynamic Linear Solver")
    sDynamic = NetPricing.ConjugateDynamicLinearModel(prob)
    tDynamic = @benchmark solve($sDynamic, d) setup=(d = $make_demands())
    display(tDynamic)
    println()

    println("KKT Solver")
    sKKT = NetPricing.ConjugateKKTModel(prob)
    tKKT = @benchmark solve($sKKT, d) setup=(d = $make_demands())
    display(tKKT)
    println()
end
