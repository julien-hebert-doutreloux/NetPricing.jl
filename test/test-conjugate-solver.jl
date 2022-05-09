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
end