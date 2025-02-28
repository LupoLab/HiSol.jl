module Focusing
import HiSol.Data: n2_solid, elastic_limit
import HiSol.Solitons: T0P0
import Luna.Maths: gauss, fwhm
import Luna: FFTW, Grid, Fields
import Luna.Modes: hquadrature
import Luna.PhysData: c
import LinearAlgebra: mul!, ldiv!
import Luna.Capillary: besselj, get_unm
import PyPlot: plt
import Printf: @sprintf

"""
    window_distance(a, λ0, energy, τfwhm, thickness; material=:SiO2, Bmax=0.2)
    window_distance(a, λ0, peakpower, thickness; material=:SiO2, Bmax=0.2)

Calculate the minimum distance between the exit of an HCF with radius `a` and a window made of `material`
with a given `thickness` such that the B-integral for a pulse at wavelength `λ0` with `energy` and duration `τfwhm`
(or a given `peakpower`) does not exceed `Bmax`.
"""
function window_distance(a, λ0, energy, τfwhm, thickness; material=:SiO2, Bmax=0.2, shape=:sech)
    _, P0 = T0P0(τfwhm, energy; shape)
    window_distance(a, λ0, P0, thickness; material, Bmax)
end

function window_distance(a, λ0, peakpower, thickness; material=:SiO2, Bmax=0.2)
    n2 = getn2(material)
    k0 = 2π/λ0
    # Bint = n2 * k0 * I0 * thickness
    Imax = Bmax/(n2*k0*thickness)
    # I0 = 2*P0/(π*proc.w0^2)
    w0min = sqrt(2peakpower/(π*Imax))
    w0HCF = 0.64a
    w0HCF >= w0min && return 0.0
    beamsize_distance(w0HCF, λ0, w0min)
end

getn2(material::Number) = material
getn2(material::Symbol) = n2_solid(material)

function window_thickness_nonlinear(a, λ0, peakpower, distance; material=:SiO2, Bmax=0.2)
    n2 = getn2(material)
    k0 = 2π/λ0
    w0window = diverged_beam(a, λ0, distance)
    Iwindow = 2peakpower/(π*w0window^2)
    Bmax/(n2*k0*Iwindow)
end

rayleigh(w0, λ) = π*w0^2/λ

"""
    mirror_distance(a, λ0, energy, LIDT; S_fluence=5)

Calculate the minimum distance between the exit of an HCF with radius `a` and a mirror
with damage threshold `LIDT` to avoid damage for a pulse at wavelength `λ0` with a given `energy`.
The maximum fluence is taken as `LIDT/S_fluence` for safety.

!!! note
    `LIDT` must be given in SI units, i.e. J/m²

# Some example values for `LIDT`:
- Newport ultrafast silver: 0.12 J/cm² = 1200 J/m²
- Thorlabs protected silver: 0.225 J/cm² = 2250 J/m²
- Layertec fs-optimised protected silver: 0.38 J/cm² = 3800 J/m²
- High-power mirror for sub-picosecond pulses: 0.75 J/cm² = 7500 J/m²
- High-LIDT mirror for Ti:Sapph pulses: 1 J/cm² = 10000 J/m²
"""
function mirror_distance(a, λ0, energy, LIDT; S_fluence=5)
    w0 = 0.64a
    zr = rayleigh(w0, λ0)
    maxfluence = LIDT/S_fluence
    sqrt(2energy*zr^2/(π*w0^2*maxfluence))
end

"""
    max_flength(a, λ0, energy, τfwhm, maxlength; kwargs...)

Find the maximum HCF length given which can fit into a space of `maxlength` given a core radius `a` and a pulse
at wavelength `λ0` with energy `energy`, FWHM duration `τfwhm`, taking into account both the damage threshold
of the end mirrors and the nonlinear lens in any windows.

The keywords `entrance_window` and `exit_window` set whether a window is present at either end.

The remaining keyword arguments
are passed to:

- [`window_distance`](@ref): `thickness`, `material`, `Bmax`
- [`mirror_distance`](@ref): `LIDT`, `S_fluence`
"""
function max_flength(a, λ0, energy, τfwhm, maxlength;
                     thickness=1e-3, material=:SiO2, Bmax=0.2,
                     LIDT=2000, S_fluence=5,
                     entrance_window=true, exit_window=true)

    dwin = window_distance(a, λ0, energy, τfwhm, thickness; material, Bmax)
    dmir = mirror_distance(a, λ0, energy, LIDT; S_fluence)

    # TODO: if both window and mirror are present, mirror is always
    # a little further away than the window
    L_in = entrance_window ? max(dwin, dmir) : dmir
    L_out = exit_window ? max(dwin, dmir) : dmir
    max(maxlength - L_in - L_out, 0)
end

"""
    diverged_beam(a, λ0, distance)

Calculate the 1/e² beam radius for a beam at wavelength `λ0` at `distance`
away from an HCF with core radius `a`.
"""
function diverged_beam(a, λ0, distance)
    w0 = 0.64a
    zr = rayleigh(w0, λ0)
    w0*sqrt(1 + (distance/zr)^2)
end

function beamsize_distance(w0, λ0, w0_diverged)
    zr = rayleigh(w0, λ0)
    zr*sqrt((w0_diverged/w0)^2 - 1)
end

function diverged_ROC(a, λ0, distance)
    w0 = 0.64a
    zr = rayleigh(w0, λ0)
    distance*(1 + (zr/distance)^2)
end

function kerr_lens(w0, energy, τfwhm, thickness; material=:SiO2, shape=:sech)
    n2 = getn2(material)
    _, P0 = T0P0(τfwhm, energy; shape)
    π*w0^4/(8*n2*thickness*P0)
end

function kerr_lens(w0, peakpower, thickness; material=:SiO2)
    n2 = getn2(material)
    π*w0^4/(8*n2*thickness*peakpower)
end

"""
    get_Bint(λ0, τfwhm, peakpower, w0, thickness; material=:SiO2, prop=true)

Calculate B-integral for a pulse at `λ0` with duration `τfwhm` and `peakpower` in a beam
with 1/e² radius `w0` when propagating through a window made of `material` of given `thickness`.

If `prop` is `true`, the dispersion of the window (and variable peak power) is taken into account
through numerical propagation.
"""
function get_Bint(λ0, τfwhm, peakpower, w0, thickness;
                  material=:SiO2, prop=true, pre_thickness=0, pre_material=:SiO2)
    n2 = getn2(material)
    k0 = 2π/λ0
    if prop
        grid = Grid.EnvGrid(thickness, λ0, (200e-9, 4e-6), 2000e-15)
        Et = sqrt.(gauss.(grid.t, fwhm=τfwhm))
        Pp_in = maximum(abs2.(Et))
        Eω = FFTW.fft(Et)
        if pre_thickness ≠ 0
            Fields.prop_material!(Eω, grid, pre_material, pre_thickness, λ0)
        end
        P0int, _ = hquadrature(0, thickness) do ti
            Eωprop = Fields.prop_material(Eω, grid, material, ti, λ0)
            Etprop = FFTW.ifft(Eωprop)
            maximum(abs2.(Etprop))/Pp_in * peakpower
        end
    else
        P0int = peakpower * thickness
    end
    I0int = 2*P0int/(π*w0^2)
    n2 * k0 * I0int
end

function plot_window_thickness_variable(a, pressure, peakpower, λ0, λmax;
                            Bmax=0.2, material=:SiO2, aperture_factor=2,
                            max_aperture_radius=10e-3)
    mindist = rayleigh(0.64a, λmax)
    maxdist = beamsize_distance(0.64a, λmax, max_aperture_radius/aperture_factor)
    distance = collect(range(mindist, maxdist, 512))

    w0win = diverged_beam.(a, λmax, distance)
    dNL = window_thickness_nonlinear.(a, λ0, peakpower, distance; material, Bmax)
    dP = window_thickness_breaking.(pressure-1, aperture_factor*w0win, material)

    idx = findfirst(eachindex(distance)) do ii
        dNL[ii] >= dP[ii]
    end

    plt.figure()
    plt.plot(distance*1e2, dNL*1e3; label="Nonlinear limit")
    plt.plot(distance*1e2, dP*1e3; label="Pressure limit")
    plt.plot(distance[idx]*1e2, dNL[idx]*1e3, "k.";
             label=@sprintf("%.2f mm thickness, %.2f cm away, %.2f mm aperture",
                            dNL[idx]*1e3, distance[idx]*1e2, 1e3aperture_factor*w0win[idx]))
    plt.xlabel("Distance (cm)")
    plt.ylabel("Window thickness (mm)")
    plt.legend()
    rax = plt.gca().twinx()
    rax.plot(distance*1e2, aperture_factor*w0win*1e3, "k--")
    rax.set_ylabel("Aperture radius (mm)")
end

function plot_window_thickness_fixed(a, pressure, peakpower, λ0, λmax, aperture_radius;
                                        Bmax=0.2, material=:SiO2, aperture_factor=2)
    mindist = rayleigh(0.64a, λmax)
    maxdist = beamsize_distance(0.64a, λmax, aperture_radius/aperture_factor)
    distance = collect(range(mindist, maxdist, 512))

    w0win = diverged_beam.(a, λmax, distance)
    dNL = window_thickness_nonlinear.(a, λ0, peakpower, distance; material, Bmax)
    dP = window_thickness_breaking(pressure-1, aperture_radius, material)

    idx = findfirst(eachindex(distance)) do ii
        dNL[ii] >= dP
    end

    plt.figure()
    plt.plot(distance*1e2, dNL*1e3; label="Nonlinear limit")
    plt.axhline(dP*1e3; color="C1", label="Pressure limit")
    plt.plot(distance[idx]*1e2, dNL[idx]*1e3, "k.";
    label=@sprintf("%.2f mm thickness, %.2f cm away",
        dNL[idx]*1e3, distance[idx]*1e2))
    plt.xlabel("Distance (cm)")
    plt.ylabel("Window thickness (mm)")
    plt.legend()
end

"""
    window_thickness_distance(a, pressure, τfwhm, energy, λ0, λmax, aperture_radius;
        Bmax=0.2, material=:SiO2, aperture_factor=2, roundmm=true, S_break=4.0,
        shape=:sech)
    window_thickness_distance(a, pressure, peakpower, λ0, λmax, aperture_radius;
        Bmax=0.2, material=:SiO2, aperture_factor=2, roundmm=true, S_break=4.0)

Calculate the required window thickness and the mimimum and maximum distance between the
exit of an HCF and the window.

Calculate the required window thickness to hold the specified `pressure` for the specified
window `material` and `aperture_radius`, and mimimum and maximum distance between the exit
of an HCF with radius `a` and the window, such that the B-integral for a pulse at
wavelength `λ0` with `energy` and duration `τfwhm` (or a given `peakpower`) does not exceed
`Bmax`, and such that the beam at wavelength `λmax` does not diverge such that the Gaussian
beam radius is less than `aperture_factor` smaller than `aperture_radius`.
If `roundmm=true`, then force the window thickness to be the next largest mm multiple. The
window thickness is calculated accounting for the safety factor `S_break`.
"""
function window_thickness_distance(a, pressure, τfwhm, energy, λ0, λmax, aperture_radius;
                                       Bmax=0.2, material=:SiO2, aperture_factor=2,
                                       roundmm=true, S_break=4.0, shape=:sech)
    _, P0 = T0P0(τfwhm, energy; shape)
    window_thickness_distance(a, pressure, P0, λ0, λmax, aperture_radius;
                                  Bmax, material, aperture_factor, roundmm, S_break)
end

function window_thickness_distance(a, pressure, peakpower, λ0, λmax, aperture_radius;
                                       Bmax=0.2, material=:SiO2, aperture_factor=2,
                                       roundmm=true, S_break=4.0)
    maximum_distance = beamsize_distance(0.64a, λmax, aperture_radius/aperture_factor)
    pressure = max(pressure - 1.0, 1.0) # make sure can handle vacuum
    thickness = window_thickness_breaking(pressure, aperture_radius, material; S_break)
    if roundmm
        thickness = ceil(thickness / 1e-3) * 1e-3
    end
    minimum_distance = window_distance(a, λ0, peakpower, thickness; material, Bmax)
    if minimum_distance > maximum_distance
        error("unable to find solution for given parameters")
    end
    (; thickness, minimum_distance, maximum_distance)
end

getmod(material::Number) = material
getmod(material::Symbol) = elastic_limit(material)

"""
    window_thickness_breaking(Δp_bar, radius, material; S_break=4.0)

Calculate the required window thickness to hold the specified `Δp_bar` pressure differential
for the specified free aperture `radius` and window `material`. The safety factor is specified
by `S_break` and defaults to 4.0. The calculation is based on an unclamped circular window.
"""
function window_thickness_breaking(Δp_bar, radius, material; S_break=4.0)
    1.06 * radius * sqrt(S_break*Δp_bar*1e5/getmod(material))
end

end