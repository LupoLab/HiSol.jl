import Test: @test, @testset
using HiSol

@testset "gauss chirps" begin
    τfwhm = 5e-15
    τfwhm2 = 30e-15

    φ2 = HiSol.Compressor.gauss_chirp(τfwhm, τfwhm2)

    @test HiSol.Compressor.gauss_chirped_duration(τfwhm, φ2) ≈ τfwhm2
end