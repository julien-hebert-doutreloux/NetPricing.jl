using BenchmarkTools
using Random

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