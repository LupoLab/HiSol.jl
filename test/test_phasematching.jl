import Test: @test, @testset
using HiSol
import Luna.Tools: pressureRDW

@testset "Finding λ vs finding pressure" begin
a = 125e-6
gas = :HeJ
τfwhm = 10e-15
energy = 200e-6
λ0 = 800e-9
λ_target = collect(range(100e-9, 400e-9, 16))

@testset "λ_target: $(λi*1e9)" for λi in λ_target
    p = HiSol.Solitons.RDW_pressure(λi, a, gas, λ0) # no nonlinear shift
    λi2 = HiSol.Solitons.RDW_wavelength(a, gas, p, λ0)
    @test isapprox(λi2, λi, rtol=1e-4)
    p = HiSol.Solitons.RDW_pressure(λi, a, gas, λ0, τfwhm, energy) # with
    λi2 = HiSol.Solitons.RDW_wavelength(a, gas, p, λ0, τfwhm, energy)
    @test isapprox(λi2, λi, rtol=1e-4)
end

end

@testset "Compare to Luna" begin
    a = 125e-6
    gas = :HeJ
    τfwhm = 10e-15
    energy = 200e-6
    λ0 = 800e-9
    λ_target = collect(range(100e-9, 400e-9, 16))

    @testset "λ_target: $(λi*1e9)" for λi in λ_target
        p = HiSol.Solitons.RDW_pressure(λi, a, gas, λ0)
        pLuna = pressureRDW(a, gas, λi, λ0)
        @test isapprox(p, pLuna, rtol=1e-3)
    end
end

@testset "RDW to ZDW conversion" begin
    a = 125e-6
    gas = :HeJ
    λ0 = 800e-9
    @testset "$pr bar" for pr in collect(range(0.2, 6, 16))
        λRDW = HiSol.Solitons.RDW_wavelength(a, gas, pr, λ0)
        λzd = HiSol.HCF.ZDW(a, gas, pr)
        λzdp = HiSol.Solitons.RDW_to_ZDW(λ0, λRDW, gas)
        @test isapprox(λzd, λzdp; rtol=1e-6)
    end

end