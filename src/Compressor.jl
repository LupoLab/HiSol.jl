module Compressor
import HiSol.Limits: critical_density, barrier_suppression_intensity
import HiSol.Solitons: T0P0
import HiSol.HCF: Aeff0, Leff, transmission, intensity_modeavg, α
import HiSol.Focusing: max_flength, diverged_beam
import HiSol.Data: n2_gas
import Luna.PhysData: pressure
import Luna.Tools: τfw_to_τ0
import Logging: @debug, @info
import Printf: @sprintf
import Roots: find_zero
import PyPlot: plt

function optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                  thickness=1e-3, material=:SiO2, Bmax=0.2,
                  LIDT=2000, S_fluence=5, S_sf=1.5, S_ion=10,
                  entrance_window=true, exit_window=true)
    factor = τfwhm_in/τfwhm_out

    ρcrit = critical_density(gas, λ0, τfwhm_in, energy; shape=:gauss)
    ρ = ρcrit/S_sf
    pr = pressure(gas, ρ)
    n2 = n2_gas(gas, pr)
    @debug(@sprintf("ρcrit: %.4e m⁻³", ρcrit))
    @debug(@sprintf("ρ: %.4e m⁻³", ρ))
    @debug(@sprintf("pressure: %.3f bar", pr))
    @debug(@sprintf("n₂: %.4e m²/W", n2))

    _, P0 = T0P0(τfwhm_in, energy; shape=:gauss)

    φm = nonlinear_phase(factor)

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
                                     thickness, material, Bmax, LIDT, S_fluence,
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

function params_maxlength(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                          thickness=1e-3, material=:SiO2, Bmax=0.2, S_sf=1.5,
                          entrance_window=true, exit_window=true, LIDT=2000, S_fluence=5)

    ρcrit = critical_density(gas, λ0, τfwhm_in, energy; shape=:gauss)
    ρ = ρcrit/S_sf
    pr = pressure(gas, ρ)
    n2 = n2_gas(gas, pr)

    _, P0 = T0P0(τfwhm_in, energy; shape=:gauss)

    broadfac_req = τfwhm_in/τfwhm_out
    φnl_req = nonlinear_phase(broadfac_req)

    k0 = 2π/λ0
    A0 = Aeff0()

    function params(a)
        maxflength = max_flength(a, λ0, energy, τfwhm_in, maxlength;
                              thickness, material, Bmax, LIDT, S_fluence,
                              entrance_window, exit_window)
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
          P0=P0, window_distance, beamsize_w0)
    end
end

function plot_optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                       thickness=1e-3, material=:SiO2, Bmax=0.2, S_sf=1.5, S_ion=10,
                       entrance_window=true, exit_window=true, LIDT=2000, S_fluence=5,
                       amin=25e-6, amax=500e-6, Na=512,
                       dot=nothing)
    factor = τfwhm_in/τfwhm_out

    f = params_maxlength(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                         thickness, material, Bmax, S_sf,
                         entrance_window, exit_window, LIDT, S_fluence)


    Isupp = barrier_suppression_intensity(gas)
    A0 = Aeff0()

    a = collect(range(amin, amax, Na))
    params = f.(a)
    broadfac = getindex.(params, :broadening_factor)
    flength = getindex.(params, :flength)
    t = getindex.(params, :transmission)
    intensity = getindex.(params, :intensity)
    window_distance = getindex.(params, :window_distance)
    Ltot = flength .+ (entrance_window ? window_distance : zero(flength)) .+ (exit_window ? window_distance : zero(flength))
    P0 = params[1].P0
    pr = params[1].pressure
    
    t[flength .== 0] .= NaN

    try
        global aopt, flopt, _, topt = optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                        thickness, material, Bmax, S_sf, S_ion, LIDT, S_fluence, entrance_window, exit_window)
        global intopt = P0/(A0*aopt^2)
        global optparams = f(aopt)
    catch e
        global aopt = missing
    end

    if ~isnothing(dot)
        dotparams = f(dot)
    end

    fig = plt.figure()
    fig.set_size_inches(15, 4)
    plt.subplot(1, 4, 1)
    plt.plot(1e6a, broadfac)
    if ~ismissing(aopt)
        plt.plot(1e6aopt, factor, "o"; color="k", label=@sprintf("%.1f μm", 1e6aopt))
    end
    if ~isnothing(dot)
        plt.plot(1e6dot, dotparams.broadening_factor, "o"; color="b", label=@sprintf("%.1f μm", 1e6dot))
    end
    plt.axhline(factor; linestyle="--", color="k", label="Required")
    plt.xlim(1e6.*extrema(a))
    plt.ylim(ymin=0)
    plt.legend(;frameon=false, fontsize=10)
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Broadening factor")

    
    plt.subplot(1, 4, 2)
    plt.plot(1e6a, flength; label="HCF")
    plt.plot(1e6a, Ltot, "--"; label="Total")
    if ~ismissing(aopt)
        Ltot_opt = (flopt
        .+ (entrance_window ? optparams.window_distance : 0)
        .+ (exit_window ? optparams.window_distance : 0)
        )
        plt.plot(1e6aopt, flopt, "o"; color="k", label=@sprintf("%.2f m", flopt))
        plt.plot(1e6aopt, Ltot_opt, "o"; fillstyle="none", color="k", label=@sprintf("%.2f m", Ltot_opt))
    end
    if ~isnothing(dot)
        Ltot_dot = (dotparams.flength
        .+ (entrance_window ? dotparams.window_distance : 0)
        .+ (exit_window ? dotparams.window_distance : 0)
        )
        plt.plot(1e6dot, dotparams.flength, "o"; color="b", label=@sprintf("%.2f m", dotparams.flength))
        plt.plot(1e6dot, Ltot_dot, "o"; fillstyle="none", color="b", label=@sprintf("%.2f m", Ltot_dot))
    end
    plt.legend(;frameon=false, fontsize=10)
    plt.xlim(1e6.*extrema(a))
    plt.ylim(ymin=0)
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Length (m)")

    plt.subplot(1, 4, 3)
    plt.plot(1e6a, 100*t)
    if ~ismissing(aopt)
        plt.plot(1e6aopt, 100*topt, "o"; color="k", label=@sprintf("%.1f %%", 100*topt))
        plt.legend(;frameon=false, fontsize=10)
    end
    if ~isnothing(dot)
        plt.plot(1e6dot, 100dotparams.transmission, "o"; color="b", label=@sprintf("%.1f %%", 100dotparams.transmission))
        plt.legend(;frameon=false, fontsize=10)
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
    end
    if ~isnothing(dot)
        intdot = P0/(A0*dot^2)
        plt.plot(1e6dot, intdot*1e-4, "o"; color="b", label=@sprintf("%.2e W/cm\$^{-2}\$", intdot*1e-4))
    end
    plt.xlim(1e6.*extrema(a))
    plt.ylim(0, 2*Isupp*1e-4/S_ion)
    plt.legend(;frameon=false, fontsize=10)
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Intensity (W/cm\$^{-2}\$)")

    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.suptitle(@sprintf("Input peak power: %.2f GW | Pressure: %.2f bar", P0*1e-9, pr))

    return fig, f
end

nonlinear_phase(broadfac) = sqrt(3*sqrt(3)/4 * (broadfac^2 - 1))
broadening_factor(φnl) = sqrt(1 + 4/3√3*φnl^2)

function gauss_chirp(τfwhm1, τfwhm2)
    τ0 = τfw_to_τ0(max(τfwhm1, τfwhm2), :gauss)
    τ0FTL = τfw_to_τ0(min(τfwhm1, τfwhm2), :gauss)
    τ0FTL*sqrt(τ0^2 - τ0FTL^2)
end

function gauss_chirped_duration(τfwhm, φ2)
    τ0 = τfw_to_τ0(τfwhm, :gauss)
    τfwhm*sqrt(φ2^2/τ0^4 + 1)
end

end