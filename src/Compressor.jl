module Compressor
import HISOL.Limits: critical_density, barrier_suppression_intensity
import HISOL.Solitons: T0P0
import HISOL.HCF: Aeff0, Leff, transmission, intensity_modeavg
import HISOL.Focusing: window_distance
import Luna.PhysData: pressure, n2_gas
import Logging: @debug, @info
import Printf: @sprintf
import Roots: find_zero

function optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                  thickness=1e-3, material=:SiO2, zr_frac=0.2, S_sf=1.5, S_ion=10)
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
        Lwin = window_distance(aguess, λ0, energy, τfwhm_in, thickness; material, zr_frac)
        global flength = maxlength - Lwin
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
    aopt = find_zero(aguess) do a
        Lwin = window_distance(a, λ0, energy, τfwhm_in, thickness; material, zr_frac)
        global flength = maxlength - Lwin
        γLeff_this = n2*k0/(A0*a^2)*Leff(flength, a, λ0)
        γLeff_this - γLeff
    end
    @debug("Optimised", 1e6*aopt)
    t = transmission(flength, aopt, λ0)
    @debug(@sprintf("Results: a = %.3f μm, %.2f m long, transmission %.4f",
                    1e6aopt, flength, t))
    return aopt, flength, t
end
end