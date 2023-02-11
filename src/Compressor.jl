module Compressor
import HiSol.Limits: critical_density, barrier_suppression_intensity
import HiSol.Solitons: T0P0
import HiSol.HCF: Aeff0, Leff, transmission, intensity_modeavg, α
import HiSol.Focusing: max_flength, diverged_beam
import HiSol.Data: n2_gas
import HiSol.GasFlow: tube_PVflow_choked, bar, slpm
import Luna.PhysData: pressure, density
import Luna.Tools: τfw_to_τ0
import Logging: @debug, @info
import Printf: @sprintf
import Roots: find_zero
import PyPlot: plt

function fix_pressure(pr, a, flength, gas, max_pressure, max_flow, ngradient)
    pr = min(pr, max_pressure)
    flow = 0.0
    if ngradient > 0
        flow = ngradient*slpm(tube_PVflow_choked(a, flength/ngradient, gas, pr*bar))
        if flow > max_flow
            pr = find_zero(pr) do prg
                flow = ngradient*slpm(tube_PVflow_choked(a, flength/ngradient,
                                                                gas, prg*bar))
                max_flow - flow
            end
        end
    end
    flow = ngradient*slpm(tube_PVflow_choked(a, flength/ngradient, gas, pr*bar))
    pr, flow
end

function optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                  thickness=1e-3, material=:SiO2, Bmax=0.2,
                  LIDT=2000, S_fluence=5, S_sf=1.5, S_ion=10,
                  entrance_window=true, exit_window=true, ngradient=0,
                  max_pressure=60.0, max_flow=1.0)
    factor = τfwhm_in/τfwhm_out

    ρcrit = critical_density(gas, λ0, τfwhm_in, energy)
    ρ = ρcrit/S_sf
    prm = pressure(gas, ρ)
    @debug(@sprintf("ρcrit: %.4e m⁻³", ρcrit))
    @debug(@sprintf("ρ: %.4e m⁻³", ρ))
    @debug(@sprintf("crit pressure: %.3f bar", pr))
    @debug(@sprintf("crit n₂: %.4e m²/W", n2))

    _, P0 = T0P0(τfwhm_in, energy)

    φm = nonlinear_phase(factor)
    φm *= (ngradient > 0 ? 3/2 : 1.0)

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
                                     thickness, material, Bmax, LIDT, S_fluence,
                                     entrance_window, exit_window)
        if flength <= 0
            @debug(@sprintf("a = %.1f μm: required window distance is too large",1e6aguess))
            aguess *= 0.99
            continue
        end
        pr, flow = fix_pressure(prm, aguess, flength, gas, max_pressure, max_flow, ngradient)
        n2 = n2_gas(gas, pr)
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
                     Decrease the compression factor or increase the maximum length,
                     maximum pressure or maximum flow.""")
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
                                     thickness, material, Bmax, LIDT, S_fluence,
                                     entrance_window, exit_window)
        global pr, flow = fix_pressure(prm, a, flength, gas, max_pressure, max_flow, ngradient)
        n2 = n2_gas(gas, pr)
        γLeff_this = n2*k0/(A0*a^2)*Leff(flength, a, λ0)
        γLeff_this - γLeff
    end
    @debug("Optimised", 1e6*aopt)
    t = transmission(flength, aopt, λ0)
    @debug(@sprintf("Results: a = %.3f μm, %.2f m long, transmission %.4f",
                    1e6aopt, flength, t))
    return aopt, flength, pr, flow, t
end

function params_maxlength(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                          thickness=1e-3, material=:SiO2, Bmax=0.2, S_sf=1.5,
                          entrance_window=true, exit_window=true, LIDT=2000, S_fluence=5,
                          ngradient=0, max_pressure=60.0, max_flow=1.0)
    _, P0 = T0P0(τfwhm_in, energy)
    broadfac_req = τfwhm_in/τfwhm_out
    φnl_req = nonlinear_phase(broadfac_req)
    φnl_req *= (ngradient > 0 ? 3/2 : 1.0)
    k0 = 2π/λ0
    A0 = Aeff0()

    function params(a)
        maxflength = max_flength(a, λ0, energy, τfwhm_in, maxlength;
                                 thickness, material, Bmax, LIDT, S_fluence,
                                 entrance_window, exit_window)
        ρcrit = critical_density(gas, λ0, τfwhm_in, energy)
        ρ = ρcrit/S_sf
        prm = pressure(gas, ρ)
        pr, flow = fix_pressure(prm, a, maxflength, gas, max_pressure, max_flow, ngradient)
        n2 = n2_gas(gas, pr)

        γthis = n2*k0/(A0*a^2)
        Leff_req = φnl_req/P0/γthis
        maxLeff = Leff(maxflength, a, λ0)
        if maxLeff > Leff_req
            α_ = α(a, λ0)
            flength = -1/α_*log(1-α_*Leff_req)
        else
            flength = maxflength
        end

        γLeff = γthis * Leff(flength, a, λ0)
        t = transmission.(flength, a, λ0)
        intensity = P0/(A0*a^2)
        φnl = P0*γLeff
        window_distance = (maxlength-maxflength)/2
        beamsize_w0 = diverged_beam(a, λ0, window_distance)
        (;φnl=P0*γLeff, transmission=t, intensity, flength,
          broadening_factor=broadening_factor(φnl), pressure=pr, n2=n2,
          P0=P0, window_distance, beamsize_w0, flow)
    end
end

function plot_optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                       thickness=1e-3, material=:SiO2, Bmax=0.2, S_sf=1.5, S_ion=10,
                       entrance_window=true, exit_window=true, LIDT=2000, S_fluence=5,
                       amin=25e-6, amax=500e-6, Na=512,
                       ngradient=0, max_pressure=60.0, max_flow=1.0)
    factor = τfwhm_in/τfwhm_out

    f = params_maxlength(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                         thickness, material, Bmax, S_sf,
                         entrance_window, exit_window, LIDT, S_fluence, ngradient,
                         max_pressure, max_flow)

    Isupp = barrier_suppression_intensity(gas)
    A0 = Aeff0()

    a = collect(range(amin, amax, Na))
    params = f.(a)
    broadfac = getindex.(params, :broadening_factor)
    flength = getindex.(params, :flength)
    t = getindex.(params, :transmission)
    intensity = getindex.(params, :intensity)
    P0 = params[1].P0
    pr = getindex.(params, :pressure)
    flow = getindex.(params, :flow)
    
    t[flength .== 0] .= NaN

    try
        global aopt, flopt, propt, flowopt, topt = optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                        thickness, material, Bmax, S_sf, S_ion, LIDT, S_fluence, entrance_window, exit_window,
                        ngradient, max_pressure, max_flow)
        global intopt = P0/(A0*aopt^2)
    catch e
        global aopt = missing
    end

    fig = plt.figure()
    fig.set_size_inches(15, 7)
    plt.subplot(2, 4, 1)
    plt.plot(1e6a, broadfac)
    if ~ismissing(aopt)
        plt.plot(1e6aopt, factor, "o"; color="k", label=@sprintf("%.1f μm", 1e6aopt))
    end
    plt.axhline(factor; linestyle="--", color="k", label="Required")
    plt.xlim(1e6.*extrema(a))
    plt.ylim(ymin=0)
    plt.legend()
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Broadening factor")

    
    plt.subplot(2, 4, 2)
    plt.plot(1e6a, flength)
    if ~ismissing(aopt)
        plt.plot(1e6aopt, flopt, "o"; color="k", label=@sprintf("%.3f m", flopt))
        plt.legend()
    end
    plt.xlim(1e6.*extrema(a))
    plt.ylim(ymin=0)
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Fibre length (m)")

    plt.subplot(2, 4, 3)
    plt.plot(1e6a, 100*t)
    if ~ismissing(aopt)
        plt.plot(1e6aopt, 100*topt, "o"; color="k", label=@sprintf("%.1f %%", 100*topt))
        plt.legend()
    end
    plt.xlim(1e6.*extrema(a))
    plt.ylim(ymin=0)
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Transmission (%)")

    plt.subplot(2, 4, 4)
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

    plt.subplot(2, 4, 5)
    plt.plot(1e6a, flow)
    plt.axhline(max_flow; linestyle="--", color="k", label="Limit")
    if ~ismissing(aopt)
        plt.plot(1e6aopt, flowopt, "o"; color="k", label=@sprintf("%.3f slpm", flowopt))
        plt.legend()
    end
    plt.xlim(1e6.*extrema(a))
    plt.ylim(ymin=0)
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Gas flow (slpm)")

    plt.subplot(2, 4, 6)
    plt.plot(1e6a, pr)
    plt.axhline(max_pressure; linestyle="--", color="k", label="Limit")
    if ~ismissing(aopt)
        plt.plot(1e6aopt, propt, "o"; color="k", label=@sprintf("%.3f bar", propt))
        plt.legend()
    end
    plt.xlim(1e6.*extrema(a))
    plt.ylim(ymin=0)
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Pressure (bar)")

    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.suptitle(@sprintf("Input peak power: %.2f GW", P0*1e-9))

    return fig, f
end

nonlinear_phase(broadfac) = sqrt(3*sqrt(3)/4 * (broadfac^2 - 1))
broadening_factor(φnl) = sqrt(1 + 4/3√3*φnl^2)

function gauss_chirp(τfwhm1, τfwhm2)
    τ0 = τfw_to_τ0(max(τfwhm1, τfwhm2), :gauss)
    τ0FTL = τfw_to_τ0(min(τfwhm1, τfwhm2), :gauss)
    τ0FTL*sqrt(τ0^2 - τ0FTL^2)
end

end