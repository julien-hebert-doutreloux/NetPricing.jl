@testset "Strong bilevel feasibility" begin
    testers = [
        ConjugateKKTTester,
        ConjugatePrimalTester,
    ]
    
    function test_same_bf(prob, paths)
        for tester_type in testers
            tester = tester_type(prob, 2)
            for p in paths, q in paths
                isbf = is_bilevel_feasible(tester, [p, q])
                issbf = is_strongly_bilevel_feasible(tester, [p, q])
                same_comm = (p[1] == q[1] && p[end] == q[end])  # Whether p and q are in the same commodity
                same_path = (p == q)                            # Whether p and q are the same path
                
                @test issbf == (same_comm ? same_path : isbf)
            end
        end
    end

    function test_strong_bf(prob, paths, nonbf, sbf)
        for tester_type in testers
            tester = tester_type(prob, 2)
            for (i, p) in enumerate(paths), (j, q) in enumerate(paths)
                isbf = is_bilevel_feasible(tester, [p, q])
                issbf = is_strongly_bilevel_feasible(tester, [p, q])
                
                @test isbf == ((i, j) ∉ nonbf && (j, i) ∉ nonbf)
                @test issbf == (i == j || (i, j) ∈ sbf || (j, i) ∈ sbf)
            end
        end
    end

    @testset "Non-overlap Network" begin
        prob = read_problem("problems/partialparallel3.json")
        paths = [
            [1,11,12,6], [1,6],
            [2,13,12,7], [2,7],
            [3,11,12,8], [3,13,12,8], [3,8],
            [4,11,12,9], [4,13,12,9], [4,9],
            [5,11,12,10], [5,13,12,10], [5,10]
        ]

        test_same_bf(prob, paths)
    end

    @testset "Overlap Partial Series" begin
        prob = read_problem("problems/bftest-partialseries.json")
        paths = [
            [1,3,5], [1,3,4,5], [1,2,3,5], [1,5],
            [3,5], [3,4,5]
        ]
        nonbf = Set([(2, 3), (3, 6)])
        sbf = Set([(1, 5), (2, 6), (3, 5), (4, 5), (4, 6)])

        test_strong_bf(prob, paths, nonbf, sbf)
    end

    @testset "Shared Parallel" begin
        prob = read_problem("problems/bftest-sharedparallel.json")
        paths = [
            [1,5,7,3], [1,5,6,7,3], [1,3],
            [2,5,7,4], [2,5,6,7,4], [2,4],
        ]
        nonbf = Set([(1, 6), (2, 6)])
        sbf = Set([(1, 4), (2, 5), (3, 4), (3, 5), (3, 6)])

        test_strong_bf(prob, paths, nonbf, sbf)
    end

    @testset "Triple commodities" begin
        prob = read_problem("problems/bftest-triplecomm.json")
        paths = [
            [1,7,8,4], [1,4],
            [2,9,10,5], [2,5],
            [3,7,8,6], [3,9,10,6], [3,6],
        ]

        test_same_bf(prob, paths)

        # Weak triple case
        for tester_type in testers
            tester2 = tester_type(prob, 2)
            tester3 = tester_type(prob, 3)
            for (i, j, k) in [(1, 4, 6), (2, 3, 5)]
                # Strongly bilevel feasible pairwise
                @test is_strongly_bilevel_feasible(tester2, paths[[i, j]])
                @test is_strongly_bilevel_feasible(tester2, paths[[j, k]])
                @test is_strongly_bilevel_feasible(tester2, paths[[i, k]])
                # Not strongly bilevel feasible overall
                @test !is_strongly_bilevel_feasible(tester3, paths[[i, j, k]])
            end
        end
    end
end