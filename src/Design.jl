module Design
import PyPlot: plt, PyCall
import Polynomials: Polynomial, roots
import Luna: PhysData
import Luna.Plotting: cmap_colours
import HiSol.Solitons: Δβwg, Δβρ, T0P0, fission_length, N, RDW_to_ZDW, τfwhm_to_T0, N_to_energy, dispersion_length, nonlinear_length, density_area_product
import HiSol.Limits: critical_intensity, barrier_suppression_intensity, Nmin, Nmax
import HiSol.HCF: intensity_modeavg, loss_length, ZDW, αbar_a, δ, fβ2, get_unm, Aeff0, dispersion, Δ
import HiSol.Focusing: max_flength
import HiSol.Data: n2_0, n2_solid

function params_maxlength(λ_target, gas, λ0, τfwhm, maxlength;
                        thickness=1e-3, material=:SiO2, Bmax=0.2,
                        LIDT=2000, S_fluence=5,
                        entrance_window=true, exit_window=true,
                        S_ion=10, S_sf=5,
                        kwargs...)
    ρasq = Δβwg(λ_target, λ0; kwargs...)/Δβρ(λ_target, gas, λ0; kwargs...)

    Icrit = critical_intensity(λ_target, gas, λ0; kwargs...)
    Isupp = barrier_suppression_intensity(gas)
    A0 = Aeff0(;kwargs...)
    n20 = n2_0(gas)
    ω0 = PhysData.wlfreq(λ0)

    zdw = RDW_to_ZDW(λ0, λ_target, gas; kwargs...)
    Nma = Nmax(zdw, gas, λ0, τfwhm; S_ion, S_sf)

    δ_ = δ(gas, λ0, zdw; kwargs...)

    function params_energy(a)
        density = ρasq/a^2
        pressure = PhysData.pressure(gas, density)
        
        Lloss = loss_length(a, λ0; kwargs...)
        
        aeff = A0*a^2
        n2 = n20 * density
        γ = ω0/PhysData.c*n2/aeff
        β2 = δ_/a^2
        
        function params(energy)
            T0, P0 = T0P0(τfwhm, energy)
            intensity = P0/aeff
            flength = max_flength(a, λ0, energy, τfwhm, maxlength;
                            thickness, material, Bmax, LIDT, S_fluence,
                            entrance_window, exit_window)
            Ld = dispersion_length(T0, β2)
            Lnl = nonlinear_length(P0, γ)
            Lfiss = sqrt(Ld*Lnl)
            Nsol = sqrt(Ld/Lnl)

            (;radius=a, density, pressure, intensity, flength, energy, τfwhm,
            N=Nsol, Nmin=Nmin(λ_target, λ0, τfwhm), Nmax=Nma,
            Lfiss, Lloss, Isupp, Icrit)
        end
    end

    params(a, energy) = params_energy(a)(energy)

    params, params_energy

end

function maxlength_limitratios(λ_target, gas, λ0, τfwhm, maxlength;
                            thickness=1e-3, material=:SiO2, Bmax=0.2,
                            entrance_window=true, exit_window=true,
                            LIDT=2000, S_fluence=5,
                            S_sf=5, S_ion=10, S_fiss=1.5, Nplot=512, kwargs...)
    _, f = params_maxlength(λ_target, gas, λ0, τfwhm, maxlength;
                            thickness, material, Bmax,
                            entrance_window, exit_window,
                            LIDT, S_fluence,
                            S_ion, S_sf, kwargs...)

    fae, f = params_maxlength(λ_target, gas, λ0, τfwhm, maxlength;
                        thickness, material, Bmax,
                        entrance_window, exit_window,
                        LIDT, S_fluence,
                        S_ion, S_sf, kwargs...)
        
    emin, amin = min_energy(λ_target, λ0, gas, τfwhm; S_sf, S_ion, S_loss=1, kwargs...)
    emax, amax = max_energy(λ_target, λ0, gas, τfwhm, maxlength;
                    S_sf, S_ion, S_fiss, thickness, material, Bmax,
                    entrance_window, exit_window, LIDT, S_fluence, kwargs...)

    a = collect(range(0.9amin, 1.5amax, Nplot))
    energy = collect(range(0.9emin, 1.1emax, Nplot))

    p = mapreduce(hcat, a) do ai
        fa = f(ai)
        map(energy) do ei
            fa(ei)
        end
    end

    Lloss = getindex.(p, :Lloss)
    flength = getindex.(p, :flength)
    Lfiss = getindex.(p, :Lfiss)

    loss_ratio = (S_fiss .* Lfiss) ./ Lloss
    fiss_ratio = (S_fiss .* Lfiss) ./flength

    N = getindex.(p, :N)
    Nmax = getindex.(p, :Nmax)
    Nmin_ratio = getindex.(p, :Nmin)./N
    Nmax_ratio = N./Nmax

    loss_idcs = (loss_ratio .< 1)
    fiss_idcs = (fiss_ratio .< 1)
    min_idcs = (Nmin_ratio .< 1)
    max_idcs = (Nmax_ratio .< 1)
    goodidcs = @. loss_idcs && fiss_idcs && min_idcs && max_idcs

    ratios = (loss=loss_ratio, fiss=fiss_ratio, Nmin=Nmin_ratio, Nmax=Nmax_ratio)
    idcs = (loss=loss_idcs, fiss=fiss_idcs, Nmin=min_idcs, Nmax=max_idcs, all=goodidcs)
    params = (;Lfiss, Lloss, N, Nmax=getindex.(p, :Nmax)[1], Nmin=getindex.(p, :Nmin)[1])

    a, energy, ratios, params, idcs, fae
end

function boundaries(
    λ_target, gas, λ0, τfwhm, maxlength;
    thickness=1e-3, material=:SiO2, Bmax=0.2,
    entrance_window=true, exit_window=true,
    LIDT=2000, S_fluence=5,
    S_sf=5, S_ion=10, S_fiss=1.5, Nplot=1024,
    kwargs...)

    emin_loss, amin_loss = min_energy(λ_target, λ0, gas, τfwhm; S_sf, S_ion, S_loss=S_fiss, kwargs...)
    emax_Nmax, amax_Nmax = max_energy(
        λ_target, λ0, gas, τfwhm, maxlength;
        thickness, material, Bmax,
        entrance_window, exit_window,
        LIDT, S_fluence,
        S_ion, S_sf, S_fiss, kwargs...)

    energy = collect(range(0.9emin_loss, 1.1emax_Nmax, Nplot))

    amax = maximum_radius.(
        λ_target, gas, λ0, τfwhm, energy, maxlength;
        thickness, material, Bmax,
        entrance_window, exit_window,
        LIDT, S_fluence,
        S_fiss, kwargs...)

    a = collect(range(0.9amin_loss, 1.5maximum(amax), Nplot))

    ρasq = density_area_product(λ_target, gas, λ0; kwargs...)
    λzd = RDW_to_ZDW(λ0, λ_target, gas; kwargs...)

    N2e = N_to_energy(gas, λ0, τfwhm; ρasq, kwargs...)

    Nmin_ = Nmin(λ_target, λ0, τfwhm)
    Nmax_ = Nmax(λzd, gas, λ0, τfwhm; S_sf, S_ion, kwargs...)

    energy_Nmin = N2e.(Nmin_, a)
    energy_Nmax = N2e.(Nmax_, a)

    ii = (energy_Nmin .< emin_loss) .&& (energy_Nmax .> emin_loss)
    bottom = hcat(a[ii], emin_loss*ones(count(ii)))

    amax_Nmax = maximum_radius.(
        λ_target, gas, λ0, τfwhm, energy_Nmax, maxlength;
        thickness, material, Bmax,
        entrance_window, exit_window,
        LIDT, S_fluence,
        S_fiss, kwargs...)
    ii = (energy_Nmax .> emin_loss) .&& (a .< amax_Nmax)
    left = hcat(a[ii], energy_Nmax[ii])

    amax_Nmin = maximum_radius.(
        λ_target, gas, λ0, τfwhm, energy_Nmin, maxlength;
        thickness, material, Bmax,
        entrance_window, exit_window,
        LIDT, S_fluence,
        S_fiss, kwargs...)
    ii = (energy_Nmin .> emin_loss) .&& (a .< amax_Nmin)
    right = hcat(a[ii], energy_Nmin[ii])

    emin_amax = N2e.(Nmin_, amax)
    emax_amax = N2e.(Nmax_, amax)
    ii = (energy .< emax_amax) .&& (energy .> emin_amax)
    top = hcat(amax[ii], energy[ii])

    full = (loss=emin_loss, Nmax=energy_Nmax, Nmin=energy_Nmin, length=amax)
    verts = vcat(bottom, right,
                 reverse(top; dims=1), reverse(left; dims=1))
    cropped = (loss=bottom, Nmax=left, Nmin=right, length=top, vertices=verts)

    a, energy, full, cropped

end

function aeplot_maxlength(λ_target, gas, λ0, τfwhm, maxlength;
                        thickness=1e-3, material=:SiO2, Bmax=0.2,
                        entrance_window=true, exit_window=true,
                        LIDT=2000, S_fluence=5,
                        S_sf=5, S_ion=10, S_fiss=1.5, Nplot=512, kwargs...)
    a, energy, ratios, params, idcs, f = maxlength_limitratios(
        λ_target, gas, λ0, τfwhm, maxlength;
        thickness, material, Bmax,
        entrance_window, exit_window,
        LIDT, S_fluence,
        S_sf, S_ion, S_fiss, Nplot, kwargs...)

    ab, energyb, full, cropped = boundaries(
        λ_target, gas, λ0, τfwhm, maxlength;
        thickness, material, Bmax,
        entrance_window, exit_window,
        LIDT, S_fluence,
        S_sf, S_ion, S_fiss, Nplot,
        kwargs...)

    patch = plt.Polygon(1e6cropped.vertices; closed=false,
                               facecolor="0.5", edgecolor="none", alpha=0.3)


    fig = plt.figure()
    fig.set_size_inches(12, 3.5)
    plt.subplot(1, 4, 1)
    plt.pcolormesh(1e6a, 1e6energy, ratios.loss; cmap="Spectral_r")
    plt.clim(0, 2)
    # plt.contour(1e6a, 1e6energy, idcs.loss, 0; colors="0.4")
    plt.contour(1e6a, 1e6energy, idcs.all, 0; colors="k")
    plt.axhline(1e6full.loss; color="0.4")
    # plt.plot(1e6allbounds[:, 1], 1e6allbounds[:, 2])
    plt.gca().add_patch(patch)
    plt.ylabel("Energy (μJ)")
    plt.xlabel("Core radius (μm)")
    plt.title("Loss")
    plt.ylim(extrema(1e6energy))
    plt.xlim(extrema(1e6a))
    plt.subplot(1, 4, 2)
    plt.pcolormesh(1e6a, 1e6energy, ratios.fiss; cmap="Spectral_r")
    plt.clim(0, 2)
    # plt.contour(1e6a, 1e6energy, idcs.fiss, 0; colors="0.4")
    plt.plot(full.length*1e6, energyb*1e6, color="0.4")
    plt.gca().set_yticklabels([])
    plt.xlabel("Core radius (μm)")
    plt.title("Fission length")
    plt.ylim(extrema(1e6energy))
    plt.xlim(extrema(1e6a))
    plt.subplot(1, 4, 3)
    plt.pcolormesh(1e6a, 1e6energy, ratios.Nmin; cmap="Spectral_r")
    plt.clim(0, 2)
    # plt.contour(1e6a, 1e6energy, idcs.Nmin, 0; colors="0.4")
    plt.plot(1e6ab, 1e6full.Nmin, color="0.4")
    plt.xlabel("Core radius (μm)")
    plt.gca().set_yticklabels([])
    plt.title("Minimum soliton order")
    plt.ylim(extrema(1e6energy))
    plt.xlim(extrema(1e6a))
    plt.subplot(1, 4, 4)
    plt.pcolormesh(1e6a, 1e6energy, ratios.Nmax; cmap="Spectral_r")
    plt.clim(0, 2)
    # plt.contour(1e6a, 1e6energy, idcs.Nmax, 0; colors="0.4")
    plt.plot(1e6ab, 1e6full.Nmax, color="0.4")
    plt.xlabel("Core radius (μm)")
    plt.gca().set_yticklabels([])
    plt.title("Maximum soliton order")
    plt.ylim(extrema(1e6energy))
    plt.xlim(extrema(1e6a))
    
    fig.tight_layout()

    Lfiss_ok = params.Lfiss[idcs.all]

    fig2 = plt.figure()
    fig2.set_size_inches(6, 3.5)
    plt.subplot(1, 2, 1)
    plt.pcolormesh(1e6a, 1e6energy, params.Lfiss)
    plt.clim(0, 1.5*maximum(Lfiss_ok))
    plt.contour(1e6a, 1e6energy, idcs.loss, 0; colors="0.4")
    plt.contour(1e6a, 1e6energy, idcs.fiss, 0; colors="0.4")
    plt.contour(1e6a, 1e6energy, idcs.all, 0; colors="k")
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Energy (μJ)")
    # plt.colorbar()
    # plt.colorbar(label="Fission length (m)")
    plt.subplot(1, 2, 2)
    plt.pcolormesh(1e6a, 1e6energy, params.N)
    plt.clim(1, 1.5maximum(params.Nmax))
    plt.contour(1e6a, 1e6energy, idcs.Nmax, 0; colors="0.4")
    plt.contour(1e6a, 1e6energy, idcs.Nmin, 0; colors="0.4")
    plt.contour(1e6a, 1e6energy, idcs.all, 0; colors="k")
    plt.xlabel("Core radius (μm)")
    # plt.colorbar(label="Soliton order")

    fig, f, a, energy, ratios
end

function aplot_energy_maxlength(λ_target, gas, λ0, τfwhm, energy, maxlength;
                                thickness=1e-3, material=:SiO2, Bmax=0.2,
                                entrance_window=true, exit_window=true,
                                LIDT=2000, S_fluence=5,
                                amin=10e-6, amax=350e-6, Na=128,
                                S_sf=5, S_ion=10, S_fiss=1.5, kwargs...)
    a = range(amin, amax, Na)
    f, _ = params_maxlength(λ_target, gas, λ0, τfwhm, maxlength;
                        thickness, material, Bmax,
                        entrance_window, exit_window,
                        LIDT, S_fluence,
                        S_ion, S_sf, kwargs...)

    p = map(a) do ai
        f(ai, energy)
    end

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
pump duration `τfwhm`, as determined by the loss limit.

Safety factors can be given as keyword arguments:
- `S_sf` (default 5): Maximum fraction of the critical power in the pump pulse
- `S_ion` (default 10): Maximum fraction of barrier-suppression intensity in the pump pulse
- `S_loss` (default 1): Maximum ratio L_loss/L_fiss where L_loss is the loss length
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

function min_energy_loss(λ_target, λ0, gas, τfwhm; S_fiss=1.5, kwargs...)
    αbar = αbar_a(λ0; kwargs...)
    ρasq = density_area_product(λ_target, gas, λ0; kwargs...)
    T0 = τfwhm_to_T0(τfwhm)
    Δ_ = Δ(gas, λ0, ρasq; kwargs...)
    
    αbar^2 * S_fiss^2 * T0^3*λ0*Aeff0(;kwargs...)/(π*n2_0(gas)*ρasq*abs(Δ_))
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
- `Bmax` (default 0.2): Maximum Kerr-lens-induced focal shift as a fraction
                           of the Rayleigh length.
- `S_fluence` (default 5): Maximum fraction of LIDT allowed on the end mirrors

Further `kwargs` `m`, `n` and `kind` determine the HCF mode.
"""
function max_energy(λ_target, λ0, gas, τfwhm, maxlength;
                    S_sf=5, S_ion=10, S_fiss=1.5,
                    thickness=1e-3, material=:SiO2, Bmax=0.2,
                    LIDT=2000, S_fluence=5,
                    entrance_window=true, exit_window=true, kwargs...)
    λzd = RDW_to_ZDW(λ0, λ_target, gas; kwargs...)
    N = Nmax(λzd, gas, λ0, τfwhm; S_sf, S_ion, kwargs...)
    T0 = τfwhm_to_T0(τfwhm)
    δ_ = δ(gas, λ0, λzd)
    f = fβ2(gas, λzd)
    n20 = n2_0(gas)
    n2w = n2_solid(material)
    u_nm = get_unm(;kwargs...)
    aeff0 = Aeff0(;kwargs...)

    # Solving S_fiss * L_fiss = L_tot - d_dwin/mir
    den1 = S_fiss*T0^2/(N*abs(δ_))
    den_win = sqrt(4*aeff0*n2w*thickness*π^3*0.64^2*f*abs(δ_)*N^2/
                (Bmax*λ0^2*n20*u_nm^2*T0^2))
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

function maximum_radius(λ_target, gas, λ0, τfwhm, energy, maxlength;
                        S_fiss=1.5,
                        thickness=1e-3, material=:SiO2, Bmax=0.2,
                        LIDT=2000, S_fluence=5,
                        entrance_window=true, exit_window=true, kwargs...)
    ρasq = density_area_product(λ_target, gas, λ0; kwargs...)
    T0, P0 = T0P0(τfwhm, energy)
    Δ_ = Δ(gas, λ0, ρasq; kwargs...)
    n20 = n2_0(gas)
    n2w = n2_solid(material)
    aeff0 = Aeff0(;kwargs...)
    maxfluence = LIDT/S_fluence
    xa = 0.64

    # solving Aa³ + Ba² + Ca + D = 0
    # with B = 0

    Cwin = sqrt(4*π^2*n2w*xa^2*thickness*P0/(λ0^3*Bmax))
    Cmir = sqrt(2energy*π*xa^2/(λ0^2*maxfluence))

    A = S_fiss*sqrt(T0^2*aeff0*λ0/(2π*P0*n20*abs(Δ_)*ρasq))
    B = 0
    C = (entrance_window ? Cwin : Cmir) + (exit_window ? Cwin : Cmir)
    D = -maxlength

    p = Polynomial([D, C, B, A])
    r = roots(p)
    real(filter(isreal, r))[1]
end


end