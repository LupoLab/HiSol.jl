module Compressor
import HiSol.Limits: critical_density, barrier_suppression_intensity
import HiSol.Solitons: T0P0
import HiSol.HCF: Aeff0, Leff, transmission, intensity_modeavg
import HiSol.Focusing: max_flength
import HiSol.Data: n2_gas
import Luna.PhysData: pressure
import Logging: @debug, @info
import Printf: @sprintf
import Roots: find_zero
import PyPlot: plt

function optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                  thickness=1e-3, material=:SiO2, zr_frac=0.2,
                  LIDT=2000, S_fluence=5, S_sf=1.5, S_ion=10,
                  entrance_window=true, exit_window=true)
    factor = τfwhm_in/τfwhm_out

    ρcrit = critical_density(gas, λ0, τfwhm_in, energy)
    ρ = ρcrit/S_sf
    pr = pressure(gas, ρ)
    n2 = n2_gas(gas, pr)
    @debug(@sprintf("ρcrit: %.4e m⁻³", ρcrit))
    @debug(@sprintf("ρ: %.4e m⁻³", ρ))
    @debug(@sprintf("pressure: %.3f bar", pr))
    @debug(@sprintf("n₂: %.4e m²/W", n2))

    _, P0 = T0P0(τfwhm_in, energy)

    # F**2 = 1 + 4/3sqrt3 φm**2
    # φm = sqrt(3sqrt3/4 *(F**2 - 1))
    φm = sqrt(3*sqrt(3)/4 * (factor^2 - 1))

    # φm = γP0Leff
    γLeff = φm/P0

    k0 = 2π/λ0
    @debug(@sprintf("k0: %.4e", k0))

    @debug(@sprintf("Required γ×Leff: %.4e", γLeff))

    Isupp = barrier_suppression_intensity(gas)
    A0 = Aeff0()
    aguess = 10e-3
    enough = false
    # Find rough guess by decreasing the core size from something massive
    # This is a reliable way of finding out whether there is a solution
    while ~enough
        global flength = max_flength(aguess, λ0, energy, τfwhm_in, maxlength;
                                     thickness, material, zr_frac, LIDT, S_fluence,
                                     entrance_window, exit_window)
        if flength <= 0
            @debug(@sprintf("a = %.1f μm: required window distance is too large",1e6aguess))
            aguess *= 0.9
            continue
        end
        γLeff_this = n2*k0/(A0*aguess^2)*Leff(flength, aguess, λ0)
        @debug(
            @sprintf(
                "a = %.1f μm: flength %.2f m | γ %.4e | γ×Leff: %.4e",
                1e6aguess, flength, n2*k0/(A0*aguess^2), γLeff_this
            )
        )

        enough = γLeff_this > γLeff
        if aguess < 10e-6
            error("""Could not find core radius: not enough nonlinearity.
                     Decrease the compression factor or increase the maximum length.""")
        end
        if P0/(A0*aguess^2) > Isupp/S_ion
            error("""Could not find core radius: intensity too high.
                     Decrease the energy or increase the maximum length.""")
        end
        aguess *= 0.9
    end
    @debug("Initial guess", 1e6*aguess)
    @debug("Initial fibre length", flength)
    aopt = find_zero(aguess) do a
        global flength = max_flength(a, λ0, energy, τfwhm_in, maxlength;
                                     thickness, material, zr_frac, LIDT, S_fluence,
                                     entrance_window, exit_window)
        γLeff_this = n2*k0/(A0*a^2)*Leff(flength, a, λ0)
        γLeff_this - γLeff
    end
    @debug("Optimised", 1e6*aopt)
    t = transmission(flength, aopt, λ0)
    @debug(@sprintf("Results: a = %.3f μm, %.2f m long, transmission %.4f",
                    1e6aopt, flength, t))
    return aopt, flength, pr, t
end

function plot_optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                       thickness=1e-3, material=:SiO2, zr_frac=0.2, S_sf=1.5, S_ion=10,
                       entrance_window=true, exit_window=true, LIDT=2000, S_fluence=5,
                       amin=10e-6, amax=500e-6, Na=512)
    factor = τfwhm_in/τfwhm_out

    ρcrit = critical_density(gas, λ0, τfwhm_in, energy)
    ρ = ρcrit/S_sf
    pr = pressure(gas, ρ)
    n2 = n2_gas(gas, pr)

    _, P0 = T0P0(τfwhm_in, energy)

    # F**2 = 1 + 4/3sqrt3 φm**2
    # φm = sqrt(3sqrt3/4 *(F**2 - 1))
    φm = sqrt(3*sqrt(3)/4 * (factor^2 - 1))

    # φm = γP0Leff
    γLeff = φm/P0

    k0 = 2π/λ0

    Isupp = barrier_suppression_intensity(gas)
    A0 = Aeff0()

    a = collect(range(amin, amax, Na))
    flength = map(a) do ai
        max_flength(ai, λ0, energy, τfwhm_in, maxlength;
                    thickness, material, zr_frac, LIDT, S_fluence,
                    entrance_window, exit_window)
    end
    γLeff_a = @. n2*k0/(A0*a^2)*Leff(flength, a, λ0)
    t = transmission.(flength, a, λ0)
    t[flength .== 0] .= NaN
    intensity = @. P0/(A0*a^2)

    try
        global aopt, flopt, _, topt = optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                        thickness, material, zr_frac, S_sf, S_ion, LIDT, S_fluence, entrance_window, exit_window)
        global intopt = P0/(A0*aopt^2)
    catch
        global aopt = missing
    end

    fig = plt.figure()
    fig.set_size_inches(15, 4)
    plt.subplot(1, 4, 1)
    plt.plot(1e6a, P0*γLeff_a)
    if ~ismissing(aopt)
        plt.plot(1e6aopt, P0*γLeff, "o"; color="k", label=@sprintf("%.1f μm", 1e6aopt))
    end
    plt.axhline(P0*γLeff; linestyle="--", color="k", label="Required")
    plt.xlim(1e6.*extrema(a))
    plt.ylim(ymin=0)
    plt.legend()
    plt.xlabel("Core radius (μm)")
    plt.ylabel("P\$_0\$ × γ × \$L_{eff}\$ (req.)")
    
    plt.subplot(1, 4, 2)
    plt.plot(1e6a, flength)
    if ~ismissing(aopt)
        plt.plot(1e6aopt, flopt, "o"; color="k", label=@sprintf("%.3f m", flopt))
        plt.legend()
    end
    plt.xlim(1e6.*extrema(a))
    plt.ylim(ymin=0)
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Fibre length (m)")

    plt.subplot(1, 4, 3)
    plt.plot(1e6a, 100*t)
    if ~ismissing(aopt)
        plt.plot(1e6aopt, 100*topt, "o"; color="k", label=@sprintf("%.1f %%", 100*topt))
        plt.legend()
    end
    plt.xlim(1e6.*extrema(a))
    plt.ylim(ymin=0)
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Transmission (%)")

    plt.subplot(1, 4, 4)
    plt.plot(1e6a, intensity*1e-4)
    plt.axhline(Isupp*1e-4/S_ion; linestyle="--", color="r", label="Limit")
    if ~ismissing(aopt)
        plt.plot(1e6aopt, intopt*1e-4, "o"; color="k", label=@sprintf("%.2e W/cm\$^{-2}\$", intopt*1e-4))
        plt.legend()
    end
    plt.xlim(1e6.*extrema(a))
    plt.ylim(0, 2*Isupp*1e-4/S_ion)
    plt.legend()
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Intensity (W/cm\$^{-2}\$)")

    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.suptitle(@sprintf("Input peak power: %.2f GW | Pressure: %.2f bar", P0*1e-9, pr))

    return fig
end
end