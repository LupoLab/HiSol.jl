module Design
import PyPlot: plt
import Luna: PhysData
import Luna.Plotting: cmap_colours
import HISOL.Solitons: Δβwg, Δβρ, T0P0, fission_length, N, RDW_to_ZDW, τfwhm_to_T0, N_to_energy
import HISOL.Limits: critical_intensity, barrier_suppression_intensity, Nmin, Nmax
import HISOL.HCF: intensity_modeavg, loss_length, ZDW, αbar_a, δ, fβ2, get_unm, Aeff0
import HISOL.Focusing: max_flength
import HISOL.Data: n2_0, n2_solid

function energy_maxlength(λ_target, gas, λ0, τfwhm, energy, maxlength;
                          thickness=1e-3, material=:SiO2, zr_frac=0.2,
                          LIDT=2000, S_fluence=5,
                          entrance_window=true, exit_window=true,
                          S_ion=10, S_sf=5,
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

        Nma = Nmax(zdw, gas, λ0, τfwhm; S_ion, S_sf)

        (;radius=a, density, pressure, intensity, flength, energy, τfwhm,
          N=Nsol, Nmin=Nmin(λ_target, λ0, τfwhm), Nmax=Nma,
          Lfiss, Lloss, Isupp, Icrit)
    end

end

function aplot_energy_maxlength(args...;
                                thickness=1e-3, material=:SiO2, zr_frac=0.2,
                                entrance_window=true, exit_window=true,
                                LIDT=2000, S_fluence = 5,
                                amin=10e-6, amax=350e-6, Na=128,
                                S_sf=5, S_ion=10, S_fiss=1.5, kwargs...)
    a = range(amin, amax, Na)
    f = energy_maxlength(args...; thickness, material, zr_frac,
                                  entrance_window, exit_window,
                                  LIDT, S_fluence,
                                  S_ion, S_sf, kwargs...)

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

"""
    min_energy(λ_target, λ0, gas, τfwhm; S_sf=5, S_ion=10, S_loss=1, kwargs...)

Find the minimum energy and core size which can be used to drive RDW emission in HCF
for an RDW target wavelength `λ_target`, pump wavelength `λ0`, fill gas `gas` and 
pump duration `τfwhm`. Safety factors can be given as keyword arguments:

- `S_sf` (default 5): Maximum fraction of the critical power in the pump pulse
- `S_ion` (default 10): Maximum fraction of barrier-suppression intensity in the pump pulse
- `S_loss` (default 1): Maximum ration L_loss/L_fiss where L_loss is the loss length
                        and L_fiss is the fission length

Further `kwargs` `m`, `n` and `kind` determine the HCF mode.
"""
function min_energy(λ_target, λ0, gas, τfwhm; S_sf=5, S_ion=10, S_loss=1, kwargs...)
    Lbar = 1/αbar_a(λ0; kwargs...)
    λzd = RDW_to_ZDW(λ0, λ_target, gas; kwargs...)
    N = Nmax(λzd, gas, λ0, τfwhm; S_sf, S_ion, kwargs...)
    δ_ = δ(gas, λ0, λzd)
    T0 = τfwhm_to_T0(τfwhm)

    a = S_loss * T0^2/(N*abs(δ_)*Lbar)
    e = N_to_energy(N, a, gas, λ0, λzd, τfwhm; kwargs...)
    e, a
end

"""
    max_energy(λ_target, λ0, gas, τfwhm, maxlength; kwargs...)

Calculate the maximum energy and core radius that can be used for an RDW
emission system for RDW wavelength `λ_target` driven by pulses at `λ0` with duration
`τfwhm` using `gas`, assuming the maximum available space is `maxlength`.

# Keyword arguments
## Entrance/exit windows
- `entrance_window`: Whether an entrance window is present (if `false`, use mirror damage threshold instead)
- `exit_window`: Same but for exit window
- `thickness` : thickness of the windows (default 1 mm)
- `material`: window material (default silica, SiO2)
- `LIDT`: laser damage threshold of the mirrors (default 2000 J/m²)

## Safety factors
- `S_sf` (default 5): Maximum fraction of the critical power in the pump pulse
- `S_ion` (default 10): Maximum fraction of barrier-suppression intensity in the pump pulse
- `S_fiss` (default 1.5): Mininum HCF length is `S_fiss` times the fission length
- `zr_frac` (default 0.2): Maximum Kerr-lens-induced focal shift as a fraction
                           of the Rayleigh length.
- `S_fluence` (default 5): Maximum fraction of LIDT allowed on the end mirrors

Further `kwargs` `m`, `n` and `kind` determine the HCF mode.
"""
function max_energy(λ_target, λ0, gas, τfwhm, maxlength;
                    S_sf=5, S_ion=10, S_fiss=1.5,
                    thickness=1e-3, material=:SiO2, zr_frac=0.2,
                    LIDT=2000, S_fluence=5,
                    entrance_window=true, exit_window=true, kwargs...)
    λzd = RDW_to_ZDW(λ0, λ_target, gas; kwargs...)
    N = Nmax(λzd, gas, λ0, τfwhm; S_sf, S_ion, kwargs...)
    T0 = τfwhm_to_T0(τfwhm)
    δ_ = δ(gas, λ0, λzd)
    f = fβ2(gas, λzd)
    n20 = n2_0(gas)
    n2w = n2_solid(material)
    u_nm = get_unm(kwargs...)
    aeff0 = Aeff0(kwargs...)

    # Solving S_fiss * L_fiss = L_tot - d_dwin/mir
    den1 = S_fiss*T0^2/(N*abs(δ_))
    den_win = sqrt(8*aeff0*n2w*thickness*π^3*0.64^2*f*abs(δ_)*N^2/
                (zr_frac*λ0^2*n20*u_nm^2*T0^2))
    den_mir = sqrt(4*aeff0*N^2*π^2*0.64^2*f*abs(δ_)/
                    (T0*λ0*LIDT/S_fluence*u_nm^2*n20))
    if entrance_window
        den2 = den_win
    else
        den2 = den_mir
    end
    if exit_window
        den2 += den_win
    else
        den2 += den_mir
    end

    a = sqrt(maxlength/(den1+den2))
    e = N_to_energy(N, a, gas, λ0, λzd, τfwhm; kwargs...)
    e, a
end


end