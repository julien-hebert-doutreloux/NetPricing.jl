using Test
using NetPricing
using Random

include("test-utils.jl")

@testset "All tests" begin
    include("test-graph.jl")
    include("test-forward-solver.jl")
    include("test-conjugate-solver.jl")
    include("test-bf.jl")
    include("test-strong-bf.jl")
end