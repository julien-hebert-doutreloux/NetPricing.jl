@testset "Shortest path fixed arcs comparison" begin
    prob = read_problem("../problems/paper/d30-01.json")
    graph = NetPricingTollBranch.build_graph(prob)
    graphmodel = NetPricingTollBranch.make_graph_model(prob)
    a1 = tolled_arcs(prob)

    function setupfunc(numfixed1)
        comm = rand(prob.K)
        fixed1 = a1[randperm(length(a1))[1:numfixed1]]
        fixed0 = setdiff(rand(a1, length(a1) ÷ 4), fixed1)
        return comm.orig, comm.dest, fixed1, fixed0
    end

    @testset "Accuracy test" begin
        for t in 1:100
            o, d, f1, f0 = setupfunc(20)
            cHungarian = NetPricingTollBranch.solve_fixedarcs_hungarian(graph, o, d, prob, f1, f0)
            cGurobi = NetPricingTollBranch.solve_fixedarcs_gurobi(graphmodel, o, d, f1, f0)
            @test cHungarian ≈ cGurobi
        end
    end

    @testset "Relaxation test" begin
        for t in 1:100
            o, d, f1, f0 = setupfunc(20)
            cNormal = NetPricingTollBranch.solve_fixedarcs_hungarian(graph, o, d, prob, f1, f0, relaxed=false)
            cRelaxed = NetPricingTollBranch.solve_fixedarcs_hungarian(graph, o, d, prob, f1, f0, relaxed=true)
            @test cRelaxed <= cNormal
        end
    end
end