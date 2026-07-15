using HiSol
using Test

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
@testset "Compressor" begin
    include("test_compressor.jl")
end
@testset "HCF" begin
    include("test_HCF.jl")
end
# must come last: loads a plotting backend, and the no-backend tests rely on
# no backend having been loaded before
@testset "Plotting" begin
    include("test_plotting.jl")
end
