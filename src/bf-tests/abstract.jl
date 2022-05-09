abstract type AbstractBFTester end

abstract type AbstractStrongBFTester <: AbstractBFTester end

# Interface
set_odpairs(tester::AbstractBFTester, commodities::AbstractVector{Commodity}) =
    set_odpairs(tester, [(comm.orig, comm.dest) for comm in commodities])

set_odpairs(tester::AbstractBFTester, ks::AbstractVector{Int}) =
    set_odpairs(tester, problem(tester).K[ks])
