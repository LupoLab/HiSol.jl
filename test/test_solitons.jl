import Test: @test, @testset
using HiSol

@testset "fission length: $(1e9λ_target) nm" for λ_target in collect(range(100e-9, 400e-9, 6))
    a = 125e-6
    gas = :HeJ
    τfwhm = 10e-15
    energy = 200e-6
    λ0 = 800e-9

    pressure = HiSol.Solitons.RDW_pressure(λ_target, a, gas, λ0)
    Lf1 = HiSol.Solitons.fission_length(a, gas, pressure, λ0, τfwhm, energy)

    N = HiSol.Solitons.N(a, gas, pressure, λ0, τfwhm, energy)
    λzd = HiSol.HCF.ZDW(a, gas, pressure)
    Lf2 = HiSol.Solitons.fission_length(a, gas, λ0, τfwhm; N, λzd)
    @test isapprox(Lf1, Lf2, rtol=1e-6)
end

@testset "soliton order <-> energy: $(1e9λ_target) nm" for λ_target in collect(range(100e-9, 400e-9, 6))
    a = 125e-6
    gas = :HeJ
    τfwhm = 10e-15
    energy = 200e-6
    λ0 = 800e-9

    pressure = HiSol.Solitons.RDW_pressure(λ_target, a, gas, λ0)
    λzd = HiSol.HCF.ZDW(a, gas, pressure)
    ρasq = HiSol.Data.density(gas, pressure)*a^2

    N = HiSol.Solitons.N(a, gas, pressure, λ0, τfwhm, energy)
    e2 = HiSol.Solitons.N_to_energy(N, a, gas, λ0, λzd, τfwhm)
    e3 = HiSol.Solitons.N_to_energy(N, a, gas, λ0, τfwhm; ρasq)

    @test isapprox(e2, energy; rtol=1e-2)
    @test isapprox(e3, energy; rtol=1e-2)
end