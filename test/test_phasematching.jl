import Test: @test, @testset
import Printf: @sprintf
using HISOL

@testset "Finding λ vs finding pressure" begin
a = 125e-6
gas = :HeJ
τfwhm = 10e-15
energy = 200e-6
λ0 = 800e-9
λ_target = collect(range(100e-9, 400e-9, 16))

@testset "λ_target: $(λi*1e9)" for λi in λ_target
    p = HISOL.Solitons.RDW_pressure(λi, a, gas, λ0) # no nonlinear shift
    λi2 = HISOL.Solitons.RDW_wavelength(a, gas, p, λ0)
    @test isapprox(λi2, λi, rtol=1e-6)
end

end