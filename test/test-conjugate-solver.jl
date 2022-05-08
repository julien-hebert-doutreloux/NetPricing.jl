@testset "Conjugate solver" begin
    prob = read_problem("problems/partialparallel3.json")
    known_objs = [
        28 19 11  8  6  6;
        18 10  6  3  3  3;
        12  5  2  1  1  1;
         7  3  0  0  0  0;
         5  2  0  0  0  0;
         5  2  0  0  0  0
    ]
    known_revs = [
         0  9 16  9  8  0;
        10 17 13 14  3  3;
        12 17 14  7  4  4;
        15 10 12  3  3  3;
         8  7  4  0  0  0;
         0  3  4  0  0  0
    ]
    paths = [
        [1,11,12,6], [1,6],
        [2,13,12,7], [2,7],
        [3,11,12,8], [3,13,12,8], [3,8],
        [4,11,12,9], [4,13,12,9], [4,9],
        [5,11,12,10], [5,13,12,10], [5,10]
    ]
    nonbf = Set([(1, 6),(1, 7),(1, 10),(1, 13),(3, 10),(3, 13),(4, 6),(5, 10),(5, 13),(6, 8),(6, 10),(6, 11),(6, 13),(8, 12),(8, 13),(9, 13)])
    solvers = [ConjugateLinearModel, ConjugateDynamicLinearModel, ConjugateKKTModel]

    @testset "Solution with priority" begin
        for solver in solvers
            cmodel = solver(prob)
            for w1 in 0:5, w2 in 0:5
                obj, t = solve(cmodel, Dict(1 => w1, 2 => w2))              # Default priority
                rev = t' * [w1, w2]
                @test obj ≈ known_objs[w1+1, w2+1] atol=1e-3
                @test rev ≈ known_revs[w1+1, w2+1] atol=1e-3                # Revenue must match
            end
        end
    end

    @testset "Solution without priority" begin
        for solver in solvers
            cmodel = solver(prob)
            for w1 in 0:5, w2 in 0:5
                obj, t = solve(cmodel, Dict(1 => w1, 2 => w2), BitSet())    # No priority
                rev = t' * [w1, w2]
                @test obj ≈ known_objs[w1+1, w2+1] atol=1e-3
                @test rev ≲ known_revs[w1+1, w2+1] atol=1e-3                # Revenue must be less
            end
        end
    end

    @testset "Bilevel feasibility test" begin
        for solver in solvers
            cmodel = solver(prob, 2)
            for (i, p) in enumerate(paths), (j, q) in enumerate(paths)
                isbf = is_bilevel_feasible(cmodel, [p, q])
                @test isbf == ((i, j) ∉ nonbf && (j, i) ∉ nonbf)
            end
        end
    end
end