using BenchmarkTools
using Random

function run_benchmark_fixedarcsolvers(probfile="problems/paper/d30-01.json")
    prob = read_problem(probfile)
    graph = NetPricingTollBranch.build_graph(prob)
    graphmodel = NetPricingTollBranch.make_graph_model(prob)
    a1 = tolled_arcs(prob)

    function setupfunc(numfixed1)
        comm = rand(prob.K)
        fixed1 = a1[randperm(length(a1))[1:min(numfixed1,end)]]
        fixed0 = setdiff(rand(a1, length(a1) ÷ 4), fixed1)
        return comm.orig, comm.dest, fixed1, fixed0
    end
    
    for nfixed1 in [1, 2, 3, 5, 7, 10, 15, 20, 30, 50]
        tHungarian = @benchmark NetPricingTollBranch.solve_fixedarcs_hungarian($graph, d[1], d[2], $prob, d[3], d[4]) setup=(d = $setupfunc($nfixed1)) evals=1
        tGurobi = @benchmark NetPricingTollBranch.solve_fixedarcs_gurobi($graphmodel, d[1], d[2], d[3], d[4]) setup=(d = $setupfunc($nfixed1)) evals=1
        @printf "Fix %2d:    %-30s " nfixed1 median(tHungarian)
        println(judge(median(tHungarian), median(tGurobi)))
    end
end

function run_benchmark_shortest_path(probfile="problems/paper/d30-01.json")
    prob = read_problem(probfile)
    graph = NetPricingTollBranch.build_graph(prob)
    
    println("Built-in Dijkstra")
    tNormal = @benchmark NetPricingTollBranch.shortest_path_old($graph, c.orig, c.dest) setup=(c = rand($prob.K))
    display(tNormal)
    println()

    println("Minimal Dijkstra")
    tMinimal = @benchmark NetPricingTollBranch.shortest_path($graph, c.orig, c.dest) setup=(c = rand($prob.K))
    display(tMinimal)
    println()

    println(judge(median(tMinimal), median(tNormal)))
end

function run_benchmark_enumeration(probfile="problems/paper/d30-01.json", numpaths=10000; kwargs...)
    prob = read_problem(probfile)
    graph = NetPricingTollBranch.build_graph(prob)
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