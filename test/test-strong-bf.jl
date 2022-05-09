@testset "Strong bilevel feasibility" begin
    prob = read_problem("problems/partialparallel3.json")
    paths = [
        [1,11,12,6], [1,6],
        [2,13,12,7], [2,7],
        [3,11,12,8], [3,13,12,8], [3,8],
        [4,11,12,9], [4,13,12,9], [4,9],
        [5,11,12,10], [5,13,12,10], [5,10]
    ]

    @testset "Pair of commodities" begin
        cmodel = ConjugateKKTModel(prob, 2)
        for p in paths, q in paths
            isbf = is_bilevel_feasible(cmodel, [p, q])
            issbf = is_strongly_bilevel_feasible(cmodel, [p, q])
            same_comm = (p[1] == q[1] && p[end] == q[end])  # Whether p and q are in the same commodity
            same_path = (p == q)                            # Whether p and q are the same path
            
            @test issbf == (same_comm ? same_path : isbf)
        end
    end
end