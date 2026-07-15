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
import HiSol: Plotting

function optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                  input_constraint, output_constraint,
                  S_sf=1.5, S_ion=10)
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
        global flength = max_flength(aguess, energy, τfwhm_in, maxlength,
                                     input_constraint, output_constraint; pressure=pr)
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
        global flength = max_flength(a, energy, τfwhm_in, maxlength,
                                     input_constraint, output_constraint; pressure=pr)
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
                          input_constraint, output_constraint, S_sf=1.5)

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
        d_in = input_constraint(a, energy, τfwhm_in; pressure=pr)
        d_out = output_constraint(a, energy, τfwhm_in; pressure=pr)
        maxflength = max(maxlength - d_in - d_out, 0)
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
        beamsize_w0_in = diverged_beam(a, λ0, d_in)
        beamsize_w0_out = diverged_beam(a, λ0, d_out)
        (;φnl=P0*γLeff, transmission=t, intensity, flength,
          broadening_factor=broadening_factor(φnl), pressure=pr, n2=n2,
          P0=P0, d_in, d_out, beamsize_w0_in, beamsize_w0_out)
    end
end

function plot_optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                       input_constraint, output_constraint,
                       S_sf=1.5, S_ion=10,
                       amin=25e-6, amax=500e-6, Na=512,
                       dot=nothing)
    factor = τfwhm_in/τfwhm_out

    f = params_maxlength(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                         input_constraint, output_constraint, S_sf)


    Isupp = barrier_suppression_intensity(gas)
    A0 = Aeff0()

    a = collect(range(amin, amax, Na))
    params = f.(a)
    broadfac = getindex.(params, :broadening_factor)
    flength = getindex.(params, :flength)
    t = getindex.(params, :transmission)
    intensity = getindex.(params, :intensity)
    d_in = getindex.(params, :d_in)
    d_out = getindex.(params, :d_out)
    Ltot = flength .+ d_in .+ d_out
    P0 = params[1].P0
    pr = params[1].pressure

    t[flength .== 0] .= NaN

    opt = try
        aopt, flopt, _, topt = optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                        input_constraint, output_constraint, S_sf, S_ion)
        optparams = f(aopt)
        (;a=aopt, flength=flopt, transmission=topt,
          intensity=P0/(A0*aopt^2),
          Ltot=flopt + optparams.d_in + optparams.d_out)
    catch e
        nothing
    end

    dotdata = if isnothing(dot)
        nothing
    else
        dotparams = f(dot)
        (;a=dot, flength=dotparams.flength, transmission=dotparams.transmission,
          intensity=P0/(A0*dot^2),
          broadening_factor=dotparams.broadening_factor,
          Ltot=dotparams.flength + dotparams.d_in + dotparams.d_out)
    end

    plotdata = (;a, broadfac, factor, flength, Ltot, transmission=t, intensity,
                 Isupp, S_ion, P0, pressure=pr,
                 opt, dot=dotdata)
    fig = Plotting.getext().plot_optimise(plotdata)

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