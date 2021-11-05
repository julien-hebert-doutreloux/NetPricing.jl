using Test
using NetPricingTollBranch
using Random

@testset "All tests" begin
    include("test-fixedarcsolver.jl")
    include("test-graph.jl")
end