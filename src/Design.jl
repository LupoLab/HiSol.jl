module Design
import PyPlot: plt
import Luna: PhysData
import Luna.PhysData: wlfreq
import Luna.Plotting: cmap_colours
import Roots: find_zero
import HiSol.Solitons: T0P0, RDW_to_ZDW, τfwhm_to_T0, N_to_energy, dispersion_length, nonlinear_length, density_area_product
import HiSol.Limits: critical_intensity, barrier_suppression_intensity, Nmin, Nmax
import HiSol.HCF: αbar_a, δ, Aeff0, Δ
import HiSol.Focusing: max_flength, LengthConstraint, NoConstraint, FixedConstraint,
                       DamageConstraint, WindowConstraint, details
import HiSol.Data: n2_0

struct Params{icT, ocT}
    λ_target::Float64
    λ0::Float64
    gas::Symbol
    maxlength::Float64
    incon::icT # length constraint on input side
    outcon::ocT # length constraint on output side
    ρasq::Float64 # density-area product
    Icrit::Float64 # critical intensity for self-focusing
    Isupp::Float64 # barrier suppression intensity
    δ::Float64 # dispersion parameter
    zdw::Float64 # zero-dispersion wavelength
    A0::Float64 # effective-area scaling
    mode::Tuple{Symbol, Int64, Int64} # mode definition: (kind, n, m)
    n20::Float64 # density-scaled nonlinear refractive index
    αbar::Float64 # core-radius-scaled attenuation coefficient
    S_ion::Float64 # safety factor on ionisation
    S_sf::Float64 # safety factor on self-focusing
    S_fiss::Float64 # safety factor on fission length
    S_loss::Float64 # safety factor on loss
end

function Params(λ_target, gas, λ0, maxlength;
                input_constraint, output_constraint,
                S_ion=10, S_sf=5, S_fiss=1.5, S_loss=1, kwargs...)
    ρasq = density_area_product(λ_target, gas, λ0; kwargs...)

    Icrit = critical_intensity(λ_target, gas, λ0; kwargs...)
    Isupp = barrier_suppression_intensity(gas)
    A0 = Aeff0(;kwargs...)
    n20 = n2_0(gas)
    αbar = αbar_a(λ0; kwargs...)
    zdw = RDW_to_ZDW(λ0, λ_target, gas; kwargs...)

    δ_ = δ(gas, λ0, zdw; kwargs...)
    kwd = Dict(kwargs)
    mode = (get(kwd, :kind, :HE), get(kwd, :n, 1), get(kwd, :m, 1))

    Params(
        float(λ_target), float(λ0), gas, float(maxlength), input_constraint, output_constraint,
        ρasq, Icrit, Isupp, δ_, zdw, A0, mode, n20, αbar,
        float(S_ion), float(S_sf), float(S_fiss), float(S_loss)
    )
end

(p::Params)(a, energy, τfwhm) = RadiusParams(a, p)(energy, τfwhm)

struct RadiusParams{pT}
    p::pT # general Params
    a::Float64 # core radius
    density::Float64 # gas density
    pressure::Float64 # gas pressure
    Lloss::Float64 # loss length
    aeff::Float64 # effective area
    n2::Float64 # nonlinear refractive index
    γ::Float64 # nonlinear coefficient
    β2::Float64 # GVD
end

function RadiusParams(a, p::Params)
    density = p.ρasq/a^2
    pressure = PhysData.pressure(p.gas, density)
    Lloss = a^3/p.αbar
    aeff = p.A0*a^2
    n2 = p.n20 * density
    γ = wlfreq(p.λ0)/PhysData.c*n2/aeff
    β2 = p.δ/a^2
    RadiusParams(p, a, density, pressure, Lloss, aeff, n2, γ, β2)
end

function (rp::RadiusParams)(energy, τfwhm)
    p = rp.p
    T0, P0 = T0P0(τfwhm, energy)
    intensity = P0/rp.aeff

    flength = max_flength(rp.a, energy, τfwhm, p.maxlength, p.incon, p.outcon;
                          pressure=rp.pressure)

    Ld = dispersion_length(T0, rp.β2)
    Lnl = nonlinear_length(P0, rp.γ)
    Lfiss = sqrt(Ld*Lnl)
    Nsol = sqrt(Ld/Lnl)
    kind, n, m = p.mode
    Nma = Nmax(p.zdw, p.gas, p.λ0, τfwhm; p.S_ion, p.S_sf, kind, n, m)
    Nmi = Nmin(p.λ_target, p.λ0, τfwhm)
    (;radius=rp.a, density=rp.density, pressure=rp.pressure, intensity, flength,
      energy, τfwhm, N=Nsol, Nmin=Nmi, Nmax=Nma, Lfiss, Lloss=rp.Lloss,
      Isupp=p.Isupp, Icrit=p.Icrit)
end

function maxlength_limitratios(λ_target, gas, λ0, τfwhm, maxlength;
                               input_constraint, output_constraint,
                               S_sf=5, S_ion=10, S_loss=1, S_fiss=1.5,
                               Nplot=512, log_e=false, log_a=false,
                               kwargs...)

    params = Params(λ_target, gas, λ0, maxlength;
                    input_constraint, output_constraint,
                    S_ion, S_sf, S_fiss, S_loss, kwargs...)

    emin, amin = min_energy(λ_target, λ0, gas, τfwhm;
                            S_sf, S_ion, S_loss, kwargs...)
    emax, amax = max_energy(λ_target, λ0, gas, τfwhm, maxlength;
                            input_constraint, output_constraint,
                            S_sf, S_ion, S_fiss, kwargs...)

    if log_e
        energy = 10 .^ collect(range(log10(0.5emin), log10(1.1emax), Nplot))
    else
        energy = collect(range(0.1emin, 1.1emax, Nplot))
    end

    if log_a
        a = 10 .^ collect(range(log10(0.9amin), log10(1.3amax), Nplot))
    else
        a = collect(range(0.9amin, 1.3amax, Nplot))
    end

    p = mapreduce(hcat, a) do ai
        rp = RadiusParams(ai, params)
        map(energy) do ei
            rp(ei, τfwhm)
        end
    end

    Lloss = getindex.(p, :Lloss)
    flength = getindex.(p, :flength)
    Lfiss = getindex.(p, :Lfiss)

    loss_ratio = (S_fiss .* Lfiss) ./ (Lloss/S_loss)
    fiss_ratio = (S_fiss .* Lfiss) ./ flength

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
    paramst = (;Lfiss, Lloss, N, Nmax=getindex.(p, :Nmax)[1], Nmin=getindex.(p, :Nmin)[1])

    a, energy, ratios, paramst, idcs, params
end

function boundaries(
    λ_target, gas, λ0, τfwhm, maxlength;
    input_constraint, output_constraint,
    S_sf=5, S_ion=10, S_fiss=1.5, Nplot=1024,
    kwargs...)

    emin_loss, amin_loss = min_energy(λ_target, λ0, gas, τfwhm; S_sf, S_ion, S_loss=S_fiss, kwargs...)
    emax_Nmax, amax_Nmax = max_energy(
        λ_target, λ0, gas, τfwhm, maxlength;
        input_constraint, output_constraint,
        S_ion, S_sf, S_fiss, kwargs...)

    energy = collect(range(0.9emin_loss, 1.1emax_Nmax, Nplot))

    amax = maximum_radius.(
        λ_target, gas, λ0, τfwhm, energy, maxlength;
        input_constraint, output_constraint,
        S_fiss, kwargs...)

    a = collect(range(0.9amin_loss, 1.5maximum(amax), Nplot))

    ρasq = density_area_product(λ_target, gas, λ0; kwargs...)
    λzd = RDW_to_ZDW(λ0, λ_target, gas; kwargs...)

    N2e = N_to_energy(gas, λ0, τfwhm; ρasq, kwargs...)

    Nmin_ = Nmin(λ_target, λ0, τfwhm)
    Nmax_ = Nmax(λzd, gas, λ0, τfwhm; S_sf, S_ion, kwargs...)

    energy_Nmin = N2e.(Nmin_, a)
    energy_Nmax = N2e.(Nmax_, a)

    amax_Nmax = maximum_radius.(
        λ_target, gas, λ0, τfwhm, energy_Nmax, maxlength;
        input_constraint, output_constraint,
        S_fiss, kwargs...)
    ii = (energy_Nmax .> emin_loss) .&& (a .< amax_Nmax)
    top = hcat(a[ii], energy_Nmax[ii])

    amax_Nmin = maximum_radius.(
        λ_target, gas, λ0, τfwhm, energy_Nmin, maxlength;
        input_constraint, output_constraint,
        S_fiss, kwargs...)
    ii = (energy_Nmin .> emin_loss) .&& (a .< amax_Nmin)
    bottom = hcat(a[ii], energy_Nmin[ii])

    emin_amax = N2e.(Nmin_, amax)
    emax_amax = N2e.(Nmax_, amax)
    ii = (energy .< emax_amax) .&& (energy .> emin_amax) .&& (energy .> emin_loss)
    right = hcat(amax[ii], energy[ii])

    amax_loss = maximum_radius.(
        λ_target, gas, λ0, τfwhm, emin_loss, maxlength;
        input_constraint, output_constraint,
        S_fiss, kwargs...)
    ii = (energy_Nmin .< emin_loss) .&& (energy_Nmax .> emin_loss) .&& (a .< amax_loss)
    left = hcat(a[ii], emin_loss*ones(count(ii)))

    full = (loss=emin_loss, Nmax=energy_Nmax, Nmin=energy_Nmin, length=amax)
    verts = vcat(bottom, right,
                 reverse(top; dims=1), reverse(left; dims=1))
    cropped = (loss=left, Nmax=top, Nmin=bottom, length=right, vertices=verts)

    a, energy, full, cropped

end

"""
    design_space_a_energy(λ_target, gas, λ0, τfwhm, maxlength;
                          input_constraint, output_constraint, kwargs...)

Map out the design space for RDW emission at `λ_target` in a `gas`-filled HCF, pumped at
wavelength `λ0` with pulses of FWHM duration `τfwhm`, where the *total* HCF system has to fit
within `maxlength`, in the plane of core radius and pulse energy.

The required keyword arguments `input_constraint` and `output_constraint` are
[`LengthConstraint`](@ref)s (e.g. [`WindowConstraint`](@ref) or [`DamageConstraint`](@ref))
which determine how much of `maxlength` is taken up by the free-space propagation on the
entrance and exit side of the HCF.

# Additional keyword arguments
- Safety factors `S_sf`, `S_ion`, `S_fiss` and `S_loss` (*larger* is always *more conservative*)
- `Nplot`: number of grid points in each direction (default: 512)
- `log_e`, `log_a`: use logarithmic spacing for the energy/radius axis (default: `false`)
- `kind`, `n`, `m`: define the HCF mode (default: HE₁₁)

Returns `(figs, params, a, energy, ratios)`, where `figs` is a tuple of `Figure`s, `params`
is a function which takes `(a, energy)` (and optionally `τfwhm` as a third argument) and
returns the full system specification at that point, `a` and `energy` are the scan axes, and
`ratios` contains the four criteria ratios on the scan grid.
"""
function design_space_a_energy(λ_target, gas, λ0, τfwhm, maxlength;
                         input_constraint, output_constraint,
                         Nplot=512, log_e=false, log_a=false, kwargs...)
    a, energy, ratios, params, idcs, f = maxlength_limitratios(
        λ_target, gas, λ0, τfwhm, maxlength;
        input_constraint, output_constraint,
        Nplot, log_e, log_a, kwargs...)

    ab, energyb, full, cropped = boundaries(
        λ_target, gas, λ0, τfwhm, maxlength;
        input_constraint, output_constraint,
        Nplot, kwargs...)

    patch = plt.Polygon(1e6cropped.vertices; closed=false,
                               facecolor="0.5", edgecolor="none", alpha=0.3)


    fig = plt.figure()
    fig.set_size_inches(12, 3.5)
    plt.subplot(1, 4, 1)
    plt.pcolormesh(1e6a, 1e6energy, ratios.loss; cmap="Spectral_r", rasterized=true)
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
    plt.pcolormesh(1e6a, 1e6energy, ratios.fiss; cmap="Spectral_r", rasterized=true)
    plt.clim(0, 2)
    # plt.contour(1e6a, 1e6energy, idcs.fiss, 0; colors="0.4")
    plt.plot(full.length*1e6, energyb*1e6, color="0.4")
    plt.gca().set_yticklabels([])
    plt.xlabel("Core radius (μm)")
    plt.title("Fission length")
    plt.ylim(extrema(1e6energy))
    plt.xlim(extrema(1e6a))
    plt.subplot(1, 4, 3)
    plt.pcolormesh(1e6a, 1e6energy, ratios.Nmin; cmap="Spectral_r", rasterized=true)
    plt.clim(0, 2)
    # plt.contour(1e6a, 1e6energy, idcs.Nmin, 0; colors="0.4")
    plt.plot(1e6ab, 1e6full.Nmin, color="0.4")
    plt.xlabel("Core radius (μm)")
    plt.gca().set_yticklabels([])
    plt.title("Minimum soliton order")
    plt.ylim(extrema(1e6energy))
    plt.xlim(extrema(1e6a))
    plt.subplot(1, 4, 4)
    plt.pcolormesh(1e6a, 1e6energy, ratios.Nmax; cmap="Spectral_r", rasterized=true)
    plt.clim(0, 2)
    # plt.contour(1e6a, 1e6energy, idcs.Nmax, 0; colors="0.4")
    plt.plot(1e6ab, 1e6full.Nmax, color="0.4")
    plt.xlabel("Core radius (μm)")
    plt.gca().set_yticklabels([])
    plt.title("Maximum soliton order")
    plt.ylim(extrema(1e6energy))
    plt.xlim(extrema(1e6a))

    fig.tight_layout()

    Lfiss_ok = copy(params.Lfiss)
    Lfiss_ok[.~idcs.all] .= NaN

    N_ok = copy(params.N)
    N_ok[.~idcs.all] .= NaN

    fig2 = plt.figure()
    fig2.set_size_inches(8, 3.5)
    gs = fig2.add_gridspec(1, 2)
    gss = gs[1].subgridspec(1, 2; width_ratios=(1, 0.05), wspace=0.05)
    ax = fig2.add_subplot(gss[1])
    img = ax.pcolormesh(1e6a, 1e6energy, Lfiss_ok, rasterized=true)
    ax.set_xlabel("Core radius (μm)")
    ax.set_ylabel("Energy (μJ)")
    ax.set_title("Fission length")
    cax = fig2.add_subplot(gss[2])
    fig2.colorbar(img; cax, label="Fission length (m)")
    gss = gs[2].subgridspec(1, 2; width_ratios=(1, 0.05), wspace=0.05)
    ax = fig2.add_subplot(gss[1])
    img = ax.pcolormesh(1e6a, 1e6energy, N_ok, rasterized=true)
    ax.set_xlabel("Core radius (μm)")
    ax.set_ylabel("Energy (μJ)")
    ax.set_title("Soliton order")
    cax = fig2.add_subplot(gss[2])
    fig2.colorbar(img; cax, label="Soliton order")
    fig2.tight_layout()

    paramsf(ai, ei, τi=τfwhm) = f(ai, ei, τi)

    (fig, fig2), paramsf, a, energy, ratios
end

function aplot_energy_maxlength(λ_target, gas, λ0, τfwhm, energy, maxlength;
                                input_constraint, output_constraint,
                                amin=10e-6, amax=350e-6, Na=128,
                                S_sf=5, S_ion=10, S_fiss=1.5, S_loss=1, kwargs...)
    a = range(amin, amax, Na)
    f = Params(λ_target, gas, λ0, maxlength;
               input_constraint, output_constraint,
               S_ion, S_sf, S_fiss, S_loss, kwargs...)

    p = map(a) do ai
        f(ai, energy, τfwhm)
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

    paramsf(ai, ei=energy, τi=τfwhm) = f(ai, ei, τi)

    return fig, paramsf
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
    δ_ = δ(gas, λ0, λzd; kwargs...)
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
    max_energy(λ_target, λ0, gas, τfwhm, maxlength; input_constraint, output_constraint, kwargs...)

Calculate the maximum energy and core radius that can be used for an RDW
emission system for RDW wavelength `λ_target` driven by pulses at `λ0` with duration
`τfwhm` using `gas`, assuming the maximum available space is `maxlength` and the free-space
propagation at either end of the HCF is determined by the [`LengthConstraint`](@ref)s
`input_constraint` and `output_constraint`.

The maximum is found at the maximum soliton order: the largest core radius (and hence energy)
for which an HCF of length `S_fiss` times the fission length, plus the distances required
by the constraints, still fits into `maxlength`.

# Keyword arguments
- `input_constraint`/`output_constraint` (required): [`LengthConstraint`](@ref)s for either end
- `S_sf` (default 5): Maximum fraction of the critical power in the pump pulse
- `S_ion` (default 10): Maximum fraction of barrier-suppression intensity in the pump pulse
- `S_fiss` (default 1.5): Minimum HCF length is `S_fiss` times the fission length

Further `kwargs` `m`, `n` and `kind` determine the HCF mode.
"""
function max_energy(λ_target, λ0, gas, τfwhm, maxlength;
                    input_constraint, output_constraint,
                    S_sf=5, S_ion=10, S_fiss=1.5,
                    kwargs...)
    λzd = RDW_to_ZDW(λ0, λ_target, gas; kwargs...)
    N = Nmax(λzd, gas, λ0, τfwhm; S_sf, S_ion, kwargs...)
    ρasq = density_area_product(λ_target, gas, λ0; kwargs...)
    N2e = N_to_energy(gas, λ0, λzd, τfwhm; kwargs...)

    a = _solve_maximum_radius(gas, λ0, τfwhm, maxlength, ρasq,
                              input_constraint, output_constraint, S_fiss;
                              kwargs...) do ai
        N2e(N, ai)
    end
    e = N2e(N, a)
    e, a
end

"""
    maximum_radius(λ_target, gas, λ0, τfwhm, energy, maxlength; input_constraint, output_constraint, kwargs...)

Calculate the maximum core radius that can be used for an RDW emission system for RDW
wavelength `λ_target` driven by pulses at `λ0` with duration `τfwhm` and `energy` using
`gas`, such that an HCF of length `S_fiss` times the fission length, plus the distances
required by `input_constraint` and `output_constraint`, still fits into `maxlength`.
Returns `0` if no core radius fits.
"""
function maximum_radius(λ_target, gas, λ0, τfwhm, energy, maxlength;
                        input_constraint, output_constraint,
                        S_fiss=1.5, kwargs...)
    ρasq = density_area_product(λ_target, gas, λ0; kwargs...)
    _solve_maximum_radius(_ -> energy, gas, λ0, τfwhm, maxlength, ρasq,
                          input_constraint, output_constraint, S_fiss;
                          err=false, kwargs...)
end

# Find the largest core radius a for which S_fiss*L_fiss + d_in + d_out == maxlength,
# where the energy at each radius is given by energyfun(a).
# If no radius fits, throw an error (err=true) or return 0.0 (err=false).
function _solve_maximum_radius(energyfun, gas, λ0, τfwhm, maxlength, ρasq,
                               input_constraint, output_constraint, S_fiss;
                               amin=1e-6, amax=20e-3, err=true, kwargs...)
    T0 = τfwhm_to_T0(τfwhm)
    Δ_ = Δ(gas, λ0, ρasq; kwargs...)
    aeff0 = Aeff0(;kwargs...)
    n20 = n2_0(gas)
    # L_fiss = prefac * a³/√P0 at constant density-area product
    Lfiss_prefac = sqrt(T0^2*aeff0*λ0/(2π*n20*abs(Δ_)*ρasq))

    function lengthdiff(a)
        energy = energyfun(a)
        _, P0 = T0P0(τfwhm, energy)
        Lf = Lfiss_prefac*a^3/sqrt(P0)
        pressure = PhysData.pressure(gas, ρasq/a^2)
        d_in = input_constraint(a, energy, τfwhm; pressure)
        d_out = output_constraint(a, energy, τfwhm; pressure)
        maxlength - S_fiss*Lf - d_in - d_out
    end

    # Coarse scan to bracket the largest radius where everything still fits. The region
    # where lengthdiff > 0 need not extend down to amin: for small cores the pressure can
    # exceed what a fixed window thickness can hold, making the constraint distance Inf.
    grid = exp10.(range(log10(amin), log10(amax), 65))
    ii = findlast(a -> lengthdiff(a) > 0, grid)
    if isnothing(ii)
        err && error(
            "Could not find a valid core radius between $(1e6amin) μm and $(1e6amax) μm: "*
            "the fission length and length constraints do not fit into "*
            "$maxlength m for any radius.")
        return 0.0
    end
    ii == lastindex(grid) && return grid[end]
    find_zero(lengthdiff, (grid[ii], grid[ii+1]))
end

end
