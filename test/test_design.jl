import Test: @test, @testset
import Roots: find_zero
import HiSol
import HiSol.Focusing: NoConstraint, FixedConstraint, DamageConstraint, WindowConstraint,
                       details, max_flength

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

# huge elastic limit disables the pressure-handling check, which the brute-force
# reference (like the old fixed-window model) does not include
window_con(λ0; thickness=1e-3, Bmax=0.2) = WindowConstraint(λ0, :SiO2;
                                                            thickness, Bmax,
                                                            elastic_limit=1e12)

@testset "max energy, $maxlength m" for maxlength in 1:20
    @testset "λ_target = $(1e9λ_target) nm" for λ_target in 100e-9:100e-9:400e-9
        λ0 = 800e-9
        gas = :HeJ
        τfwhm = 10e-15

        @testset "Bmax = $Bmax" for Bmax in 0.1:0.1:1
            wc = window_con(λ0; Bmax)
            emax, amax = HiSol.Design.max_energy(λ_target, λ0, gas, τfwhm, maxlength;
                                                 input_constraint=wc, output_constraint=wc)
            amax2 = max_energy_brute_window(λ_target, λ0, gas, τfwhm, maxlength; Bmax)
            @test isapprox(amax, amax2; rtol=1e-3)
        end
        @testset "LIDT = $LIDT" for LIDT in 1000:500:3000
            dc = DamageConstraint(λ0, LIDT)
            emax, amax = HiSol.Design.max_energy(λ_target, λ0, gas, τfwhm, maxlength;
                                                 input_constraint=dc, output_constraint=dc)
            amax2 = max_energy_brute_mirror(λ_target, λ0, gas, τfwhm, maxlength; LIDT)
            @test isapprox(amax, amax2; rtol=1e-3)
        end
    end
end

@testset "max radius, $(1e6energy) μJ, $(λ_target*1e9) nm" for energy in 1e-6*(100:100:2000), λ_target in 100e-9:100e-9:400e-9
    λ0 = 800e-9
    gas = :HeJ
    τfwhm = 10e-15
    maxlength = 5

    LIDT = 2000
    S_fluence = 5
    Bmax = 0.2
    material = :SiO2
    thickness = 1e-3
    S_fiss = 1.5

    wc = window_con(λ0; thickness, Bmax)
    dc = DamageConstraint(λ0, LIDT; S_fluence)

    @testset "input: $in_type, output: $out_type" for in_type in (:window, :mirror), out_type in (:window, :mirror)
        input_constraint = in_type == :window ? wc : dc
        output_constraint = out_type == :window ? wc : dc
        amax = HiSol.Design.maximum_radius(λ_target, gas, λ0, τfwhm, energy, maxlength;
                                           input_constraint, output_constraint, S_fiss)

        dwin = HiSol.Focusing.window_distance(amax, λ0, energy, τfwhm, thickness; material, Bmax)
        dmir = HiSol.Focusing.mirror_distance(amax, λ0, energy, LIDT; S_fluence)

        pressure = HiSol.Solitons.RDW_pressure(λ_target, amax, gas, λ0)
        Lf = HiSol.Solitons.fission_length(amax, gas, pressure, λ0, τfwhm, energy)

        d = (in_type == :window ? dwin : dmir) + (out_type == :window ? dwin : dmir)

        @test isapprox(S_fiss*Lf + d, maxlength, rtol=1e-3)
    end
end

@testset "Length constraints" begin
    a = 125e-6
    λ0 = 800e-9
    energy = 200e-6
    τfwhm = 10e-15
    pressure = 2.0

    @testset "NoConstraint and FixedConstraint" begin
        @test NoConstraint()(a, energy, τfwhm; pressure) == 0
        @test details(NoConstraint(), a, energy, τfwhm; pressure).distance == 0
        @test FixedConstraint(0.3)(a, energy, τfwhm; pressure) == 0.3
        @test FixedConstraint(1)(a, energy, τfwhm) === 1.0
        @test details(FixedConstraint(0.3), a, energy, τfwhm).distance == 0.3
    end

    @testset "DamageConstraint" begin
        dc = DamageConstraint(λ0, 2000)
        @test dc(a, energy, τfwhm) == HiSol.Focusing.mirror_distance(a, λ0, energy, 2000; S_fluence=5)
        dc2 = DamageConstraint(λ0, 2000; S_fluence=2, conversion=0.5)
        @test dc2(a, energy, τfwhm) == HiSol.Focusing.mirror_distance(a, λ0, 0.5energy, 2000; S_fluence=2)
        det = details(dc, a, energy, τfwhm)
        @test det.distance == dc(a, energy, τfwhm)
        @test det.beam_radius == HiSol.Focusing.diverged_beam(a, λ0, det.distance)
    end

    @testset "WindowConstraint: fixed thickness" begin
        wc = window_con(λ0)
        @test wc(a, energy, τfwhm; pressure) == HiSol.Focusing.window_distance(
            a, λ0, energy, τfwhm, 1e-3; material=:SiO2, Bmax=0.2)
        det = details(wc, a, energy, τfwhm; pressure)
        @test det.distance == wc(a, energy, τfwhm; pressure)
        @test det.thickness == 1e-3
        @test det.aperture == 2*HiSol.Focusing.diverged_beam(a, λ0, det.distance)
        # tiny elastic limit: window cannot hold the pressure -> Inf
        wc_weak = WindowConstraint(λ0, :SiO2; thickness=1e-3, elastic_limit=1.0)
        @test wc_weak(a, energy, τfwhm; pressure) == Inf
    end

    @testset "WindowConstraint: variable thickness and aperture" begin
        @testset "round_thickness: $rt, round_aperture: $ra" for rt in (false, true, 0.5), ra in (false, true, 0.5)
            wc = WindowConstraint(λ0, :SiO2; round_thickness=rt, round_aperture=ra)
            # window_distance_thickness_aperture self-verifies and errors if inconsistent
            det = details(wc, a, energy, τfwhm; pressure)
            @test det.distance == wc(a, energy, τfwhm; pressure)
            @test all(isfinite, det)
            @test det.distance >= 0
            @test det.thickness > 0
            @test det.aperture > 0
            rt == true && @test isinteger(round(det.thickness/1e-3; digits=6))
            rt == 0.5 && @test isinteger(round(det.thickness/0.5e-3; digits=6))
            ra == true && @test isinteger(round(det.aperture/1e-3; digits=6))
            ra == 0.5 && @test isinteger(round(det.aperture/0.5e-3; digits=6))
        end
    end

    @testset "WindowConstraint: fixed aperture" begin
        wc = WindowConstraint(λ0, :SiO2; aperture=10e-3)
        det = details(wc, a, energy, τfwhm; pressure)
        @test det.aperture == 10e-3
        @test det.thickness == HiSol.Focusing.window_thickness_breaking(
            max(pressure-1, 1), 10e-3, :SiO2; S_break=4)
        @test det.distance == HiSol.Focusing.window_distance(
            a, λ0, energy, τfwhm, det.thickness; material=:SiO2, Bmax=0.2)
        # tiny aperture: beam does not fit through -> Inf
        wc_small = WindowConstraint(λ0, :SiO2; aperture=0.1e-3)
        @test wc_small(a, energy, τfwhm; pressure) == Inf
    end

    @testset "max_flength" begin
        maxlength = 5
        fc = FixedConstraint(0.5)
        dc = DamageConstraint(λ0, 2000)
        @test max_flength(a, energy, τfwhm, maxlength, fc, fc) == 4.0
        @test max_flength(a, energy, τfwhm, maxlength, fc, dc; pressure) ==
            maxlength - 0.5 - dc(a, energy, τfwhm)
        # constraints which do not fit clamp to zero
        @test max_flength(a, energy, τfwhm, maxlength, FixedConstraint(3), FixedConstraint(3)) == 0
    end
end

@testset "Params" begin
    λ_target = 160e-9
    gas = :He
    λ0 = 800e-9
    τfwhm = 10e-15
    maxlength = 5
    a = 125e-6
    energy = 200e-6

    wc = window_con(λ0)
    p = HiSol.Design.Params(λ_target, gas, λ0, maxlength;
                            input_constraint=wc, output_constraint=wc)
    r = p(a, energy, τfwhm)
    for k in (:pressure, :density, :intensity, :N, :Nmin, :Nmax, :Lfiss, :Lloss, :flength)
        @test isfinite(getproperty(r, k))
        @test getproperty(r, k) > 0
    end
    @test r.flength ≈ maxlength - 2*wc(a, energy, τfwhm; pressure=r.pressure)

    # a RadiusParams object gives the same results as the Params object
    rp = HiSol.Design.RadiusParams(a, p)
    @test rp(energy, τfwhm) == r

    # flength clamps at zero when the constraints do not fit
    pfix = HiSol.Design.Params(λ_target, gas, λ0, maxlength;
                               input_constraint=FixedConstraint(3),
                               output_constraint=FixedConstraint(3))
    @test pfix(a, energy, τfwhm).flength == 0

    # safety factors are passed through to the soliton-order limits
    p2 = HiSol.Design.Params(λ_target, gas, λ0, maxlength;
                             input_constraint=wc, output_constraint=wc,
                             S_ion=40, S_sf=20)
    @test p2(a, energy, τfwhm).Nmax < r.Nmax

    # mode keywords are accepted and change the result
    pm = HiSol.Design.Params(λ_target, gas, λ0, maxlength;
                             input_constraint=wc, output_constraint=wc, m=2)
    @test pm.mode == (:HE, 1, 2)
    @test pm(a, energy, τfwhm).Nmax != r.Nmax
end
