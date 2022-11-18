import Test: @test, @testset
import Roots: find_zero
import HiSol

function max_energy_brute_window(λ_target, λ0, gas, τfwhm, maxlength;
                                 S_sf=5, S_ion=10, S_fiss=1.5,
                                 thickness=1e-3, material=:SiO2, Bmax=0.2,
                                 kwargs...)
    λzd = HiSol.Solitons.RDW_to_ZDW(λ0, λ_target, gas; kwargs...)
    N = HiSol.Limits.Nmax(λzd, gas, λ0, τfwhm; S_sf, S_ion, kwargs...)

    amax = find_zero((10e-6, 2e-3)) do a
        energy = HiSol.Solitons.N_to_energy(N, a, gas, λ0, λzd, τfwhm; kwargs...)
        dwin = HiSol.Focusing.window_distance(a, λ0, energy, τfwhm, thickness; material, Bmax)
        pressure = HiSol.Solitons.RDW_pressure(λ_target, a, gas, λ0; kwargs...)
        Lf = HiSol.Solitons.fission_length(a, gas, pressure, λ0, τfwhm, energy)
        LHCF = S_fiss*Lf
        maxlength - LHCF - 2*dwin
    end
    return amax
end

function max_energy_brute_mirror(λ_target, λ0, gas, τfwhm, maxlength;
                                 S_sf=5, S_ion=10, S_fiss=1.5,
                                 LIDT=2000, S_fluence=5,
                                 kwargs...)
    λzd = HiSol.Solitons.RDW_to_ZDW(λ0, λ_target, gas; kwargs...)
    N = HiSol.Limits.Nmax(λzd, gas, λ0, τfwhm; S_sf, S_ion, kwargs...)

    amax = find_zero((10e-6, 2e-3)) do a
        energy = HiSol.Solitons.N_to_energy(N, a, gas, λ0, λzd, τfwhm; kwargs...)
        dmir = HiSol.Focusing.mirror_distance(a, λ0, energy, LIDT; S_fluence)
        pressure = HiSol.Solitons.RDW_pressure(λ_target, a, gas, λ0; kwargs...)
        Lf = HiSol.Solitons.fission_length(a, gas, pressure, λ0, τfwhm, energy)
        LHCF = S_fiss*Lf
        maxlength - LHCF - 2*dmir
    end
    return amax
end

@testset "max energy, $maxlength m" for maxlength in 1:20
    @testset "λ_target = $(1e9λ_target) nm" for λ_target in 100e-9:100e-9:400e-9
        λ0 = 800e-9
        gas = :HeJ
        τfwhm = 10e-15

        @testset "Bmax = $Bmax" for Bmax in 0.1:0.1:1
            emax, amax = HiSol.Design.max_energy(λ_target, λ0, gas, τfwhm, maxlength; Bmax)
            amax2 = max_energy_brute_window(λ_target, λ0, gas, τfwhm, maxlength; Bmax)
            @test isapprox(amax, amax2; rtol=1e-3)
        end
        @testset "LIDT = $LIDT" for LIDT in 1000:500:3000
            emax, amax = HiSol.Design.max_energy(λ_target, λ0, gas, τfwhm, maxlength; LIDT, entrance_window=false, exit_window=false)
            amax2 = max_energy_brute_mirror(λ_target, λ0, gas, τfwhm, maxlength; LIDT)
            @test isapprox(amax, amax2; rtol=1e-3)
        end
    end
end