@testset "Bilevel feasibility" begin
    prob = read_problem("problems/partialparallel3.json")
    paths = [
        [1,11,12,6], [1,6],
        [2,13,12,7], [2,7],
        [3,11,12,8], [3,13,12,8], [3,8],
        [4,11,12,9], [4,13,12,9], [4,9],
        [5,11,12,10], [5,13,12,10], [5,10]
    ]
    nonbf = Set([(1, 6),(1, 7),(1, 10),(1, 13),(3, 10),(3, 13),(4, 6),(5, 10),(5, 13),(6, 8),(6, 10),(6, 11),(6, 13),(8, 12),(8, 13),(9, 13)])

    testers = [
        ConjugateBFTester{ConjugateLinearModel},
        ConjugateBFTester{ConjugateDynamicLinearModel},
        ConjugateBFTester{ConjugateKKTModel},
        ConjugateKKTTester,
        ConjugatePrimalTester,
    ]

    @testset "Pair of paths" begin
        for tester_type in testers
            tester = tester_type(prob, 2)
            for (i, p) in enumerate(paths), (j, q) in enumerate(paths)
                isbf = is_bilevel_feasible(tester, [p, q])
                @test isbf == ((i, j) ∉ nonbf && (j, i) ∉ nonbf)
            end
        end
    end
end
