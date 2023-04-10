@testset "Forward solver" begin
    prob = read_problem("problems/partialparallel3.json")
    pprobs = preprocess(prob)
    known_objs = [
        0  2  4  5  5  5  5  5  5  5  5;
        3  5  7  9  9  9  9  9  9  9  9;
        5  8 10 12 13 13 13 13 13 13 13;
        6  9 12 14 15 16 16 16 16 16 16;
        6 10 13 16 17 18 19 19 19 19 19;
        6 10 14 17 19 20 21 22 22 22 22;
        6 10 14 17 19 21 22 23 24 24 24;
        6 10 14 17 19 21 23 24 25 25 25;
        6 10 14 17 19 21 23 25 26 26 26;
        6 10 14 17 19 21 23 25 27 27 27;
        6 10 14 17 19 21 23 25 27 28 28
    ]
    known_revs = [
        0  2  4  3  0  0  0  0  0  0  0;
        3  5  7  9  4  4  4  4  4  4  4;
        4  8 10 12 10  8  8  8  8  8  8;
        3  6 10 12 10 11  9  9  9  9  9;
        0  7 10 14 12 13 14 12 12 12 12;
        0  4 11 14 14 15 16 17 15 15 15;
        0  4  8  9  8 11 12 13 14 12 12;
        0  4  8  9  8 10 13 14 15  7  7;
        0  4  8  9  8 10 12 15 16  8  8;
        0  4  8  9  8 10 12 14 17  9  9;
        0  4  8  9  8 10 12 14 16 10 10
    ]
    solvers = [
        () -> ForwardHybridSolver(pprobs; maxpaths=0),
        () -> ForwardHybridSolver(pprobs)]

    @testset "Solution with priority" begin
        for solver in solvers
            fmodel = solver()
            for t1 in 0:10, t2 in 0:10
                obj, w = solve(fmodel, Dict(1 => t1, 2 => t2))              # Default priority
                rev = w' * [t1, t2]
                @test obj ≈ known_objs[t1+1, t2+1] atol=1e-3
                @test rev ≈ known_revs[t1+1, t2+1] atol=1e-3                # Revenue must match
            end
        end
    end

    @testset "Solution without priority" begin
        for solver in solvers
            fmodel = solver()
            for t1 in 0:10, t2 in 0:10
                obj, w = solve(fmodel, Dict(1 => t1, 2 => t2), BitSet())    # No priority
                rev = w' * [t1, t2]
                @test obj ≈ known_objs[t1+1, t2+1] atol=1e-3
                @test rev ≲ known_revs[t1+1, t2+1] atol=1e-3                # Revenue must be less
            end
        end
    end
end