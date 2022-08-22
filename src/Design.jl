module Design
import PyPlot: plt
import Luna: PhysData
import Luna.Plotting: cmap_colours
import HISOL.Solitons: Δβwg, Δβρ, T0P0, fission_length, N, RDW_to_ZDW, τfwhm_to_T0, N_to_energy
import HISOL.Limits: critical_intensity, barrier_suppression_intensity, Nmin, Nmax
import HISOL.HCF: intensity_modeavg, loss_length, ZDW, αbar_a, δ
import HISOL.Focusing: max_flength

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

    loss_idcs = (loss_ratio .> 1)
    fiss_idcs = (fiss_ratio .> 1)
    min_idcs = (Nmin_ratio .> 1)
    max_idcs = (Nmax_ratio .> 1)
    goodidcs = @. loss_idcs & fiss_idcs & min_idcs & max_idcs
    agood = a[goodidcs]

    cols = cmap_colours(4)
    fig = plt.figure()
    l, = plt.plot(1e6a, loss_ratio; label="\$L_l / 1.5L_f\$", c=cols[1])
    plt.fill_between(1e6a[loss_idcs], ones(length(a[loss_idcs])), loss_ratio[loss_idcs]; color=l.get_color(), alpha=0.25)
    l, = plt.plot(1e6a, fiss_ratio; label="\$L_{max} / 1.5L_f\$", c=cols[2])
    plt.fill_between(1e6a[fiss_idcs], ones(length(a[fiss_idcs])), fiss_ratio[fiss_idcs]; color=l.get_color(), alpha=0.25)
    l, = plt.plot(1e6a, Nmin_ratio; label="\$N / N_{min}\$", c=cols[3])
    plt.fill_between(1e6a[min_idcs], ones(length(a[min_idcs])), Nmin_ratio[min_idcs]; color=l.get_color(), alpha=0.25)
    l, = plt.plot(1e6a, Nmax_ratio; label="\$N_{max} / N\$", c=cols[4])
    plt.fill_between(1e6a[max_idcs], ones(length(a[max_idcs])), Nmax_ratio[max_idcs]; color=l.get_color(), alpha=0.25)
    plt.axhline(1; color="0.5")
    if length(agood) > 1
        plt.fill_between(1e6*agood, 0, 1, color="r", alpha=0.2)
    end
    plt.legend()
    plt.ylim(0, 1.5*loss_ratio[1])
    plt.xlim(1e6.*extrema(a))
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Ratio")

    return fig, f
end

function min_energy(λ_target, λ0, gas, τfwhm; S_sf=5, S_ion=10, S_fiss=1, S_loss=1, kwargs...)
    Lbar = 1/αbar_a(λ0; kwargs...)
    λzd = RDW_to_ZDW(λ0, λ_target, gas; kwargs...)
    N = Nmax(λzd, gas, λ0, τfwhm; S_sf, S_ion, kwargs...)
    δ_ = δ(gas, λ0, λzd)
    T0 = τfwhm_to_T0(τfwhm)

    a = S_loss * S_fiss * T0^2/(N*abs(δ_)*Lbar)
    e = N_to_energy(N, a, gas, λ0, λzd, τfwhm; kwargs...)
    e, a
end


end