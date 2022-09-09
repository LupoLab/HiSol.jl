using HiSol
using Test

@testset "HiSol.jl" begin
    @testset "Phase-matching" begin
        include("test_phasematching.jl")
    end
    @testset "Limits" begin
        include("test_limits.jl")
    end
    @testset "Design" begin
        include("test_design.jl")
    end
    @testset "Solitons" begin
        include("test_solitons.jl")
    end
end
