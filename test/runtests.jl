using Test
using NetPricingTollBranch
using BenchmarkTools
using Random

@testset "Shortest path fixed arcs comparison" begin
    prob = read_problem("../problems/d30-01.json")
    graph = NetPricingTollBranch.build_graph(prob)
    graphmodel = NetPricingTollBranch.make_graph_model(prob)
    a1 = tolled_arcs(prob)

    function setupfunc(numfixed1)
        comm = rand(prob.K)
        fixed1 = a1[randperm(length(a1))[1:numfixed1]]
        fixed0 = setdiff(rand(a1, length(a1) ÷ 4), fixed1)
        return comm.orig, comm.dest, fixed1, fixed0
    end

    @testset "Accuracy tests" begin
        for t in 1:100
            o, d, f1, f0 = setupfunc(20)
            cHungarian = NetPricingTollBranch.solve_fixedarcs_hungarian(graph, o, d, prob, f1, f0)
            cGurobi = NetPricingTollBranch.solve_fixedarcs_gurobi(graphmodel, o, d, f1, f0)
            @test cHungarian ≈ cGurobi
        end
    end
    
    @testset "Benchmark" begin
        for nfixed1 in [1, 2, 3, 5, 7, 10, 15, 20, 30, 50]
            tHungarian = @benchmark NetPricingTollBranch.solve_fixedarcs_hungarian($graph, d[1], d[2], $prob, d[3], d[4]) setup=(d = $setupfunc($nfixed1)) evals=1
            tGurobi = @benchmark NetPricingTollBranch.solve_fixedarcs_gurobi($graphmodel, d[1], d[2], d[3], d[4]) setup=(d = $setupfunc($nfixed1)) evals=1
            println("Fix $nfixed1:    ", median(tHungarian), "    ", judge(median(tHungarian), median(tGurobi)))
        end
    end
end