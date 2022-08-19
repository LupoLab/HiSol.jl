module Design
import PyPlot: plt
import Luna: PhysData
import HISOL.Solitons: Δβwg, Δβρ, T0P0, fission_length, N, RDW_to_ZDW
import HISOL.Limits: critical_intensity, barrier_suppression_intensity, Nmin, Nmax
import HISOL.HCF: intensity_modeavg, loss_length, ZDW
import HISOL.Focusing: window_distance, mirror_distance

function energy_maxlength(λ_target, gas, λ0, τfwhm, energy, maxlength;
                          thickness=1e-3, material=:SiO2, zr_frac=0.2,
                          LIDT=2000, S_fluence=5,
                          entrance_window=true, exit_window=true,
                          kwargs...)
    ρasq = Δβwg(λ_target, λ0; kwargs...)/Δβρ(λ_target, gas, λ0; kwargs...)

    Icrit = critical_intensity(λ_target, gas, λ0; kwargs...)
    Isupp = barrier_suppression_intensity(gas)

    _, P0 = T0P0(τfwhm, energy)

    function params(a)
        density = ρasq/a^2
        pressure = PhysData.pressure(gas, density)
        intensity = intensity_modeavg(a, P0; kwargs...)

        global flength = max_flength(a, λ0, energy, τfwhm, maxlength;
                                     thickness, material, zr_frac, LIDT, S_fluence,
                                     entrance_window, exit_window)
        Nsol = N(a, gas, pressure, λ0, τfwhm, energy; kwargs...)

        Lfiss = fission_length(a, gas, pressure, λ0, τfwhm, energy; kwargs...)
        Lloss = loss_length(a, λ0; kwargs...)

        zdw = RDW_to_ZDW(λ0, λ_target, gas; kwargs...)

        (;radius=a, density, pressure, intensity, flength, energy, τfwhm,
          N=Nsol, Nmin=Nmin(λ_target, λ0, τfwhm), Nmax=Nmax(zdw, gas, λ0, τfwhm),
          Lfiss, Lloss, Isupp, Icrit)
    end

end

function aplot_energy_maxlength(args...;
                                exit_window=true, thickness=1e-3, material=:SiO2, zr_frac=0.2,
                                amin=10e-6, amax=350e-6, Na=128,
                                S_sf=5, S_ion=10, S_fiss=1.5, kwargs...)
    a = range(amin, amax, Na)
    f = energy_maxlength(args...; exit_window, thickness, material, zr_frac, kwargs...)

    p = map(f, a)

    Lloss = getindex.(p, :Lloss)
    flength = getindex.(p, :flength)
    Lfiss = getindex.(p, :Lfiss)

    loss_ratio = Lloss ./ (S_fiss .* Lfiss)
    fiss_ratio = flength ./ (S_fiss .* Lfiss)

    Nmin_ratio = getindex.(p, :N)./getindex.(p, :Nmin)
    Nmax_ratio = getindex.(p, :Nmax)./getindex.(p, :N)

    goodidcs = @. (loss_ratio > 1) & (fiss_ratio > 1) & (Nmin_ratio > 1) & (Nmax_ratio > 1)
    agood = a[goodidcs]

    fig = plt.figure()
    plt.plot(1e6a, loss_ratio; label="\$L_l / 1.5L_f\$")
    plt.plot(1e6a, fiss_ratio; label="\$L_{max} / 1.5L_f\$")
    plt.plot(1e6a, Nmin_ratio; label="\$N / N_{min}\$")
    plt.plot(1e6a, Nmax_ratio; label="\$N_{max} / N\$")
    plt.axhline(1; color="0.5")
    if length(agood) > 1
        plt.fill_between(1e6*agood, 0, 1, color="r", alpha=0.2)
    end
    plt.legend()
    plt.ylim(0, 5)
    plt.xlim(1e6.*extrema(a))
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Ratio")

    return fig, f
end
end