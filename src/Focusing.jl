module Focusing
import HiSol.Data: n2_solid
import HiSol.Solitons: T0P0
import Luna.Maths: gauss, fwhm
import Luna: FFTW, Grid, Fields
import Luna.Modes: hquadrature

"""
    window_distance(a, λ0, energy, τfwhm, thickness; material=:SiO2, Bmax=0.2)
    window_distance(a, λ0, peakpower, thickness; material=:SiO2, Bmax=0.2)

Calculate the minimum distance between the exit of an HCF with radius `a` and a window made of `material`
with a given `thickness` such that the B-integral for a pulse at wavelength `λ0` with `energy` and duration `τfwhm`
(or a given `peakpower`) does not exceed `Bmax`.
"""
function window_distance(a, λ0, energy, τfwhm, thickness; material=:SiO2, Bmax=0.2)
    _, P0 = T0P0(τfwhm, energy)
    window_distance(a, λ0, P0, thickness; material, Bmax)
end

function window_distance(a, λ0, peakpower, thickness; material=:SiO2, Bmax=0.2)
    n2 = n2_solid(material)
    k0 = 2π/λ0
    # Bint = n2 * k0 * I0 * thickness
    Imax = Bmax/(n2*k0*thickness)
    # I0 = 2*P0/(π*proc.w0^2)
    w0min = sqrt(2peakpower/(π*Imax))
    w0HCF = 0.64a
    w0HCF >= w0min && return 0.0
    beamsize_distance(w0HCF, λ0, w0min)
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

function kerr_lens(w0, energy, τfwhm, thickness; material=:SiO2)
    n2 = n2_solid(material)
    _, P0 = T0P0(τfwhm, energy)
    π*w0^4/(8*n2*thickness*P0)
end

function kerr_lens(w0, peakpower, thickness; material=:SiO2)
    n2 = n2_solid(material)
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
                  material=:SiO2, prop=true)
    n2 = n2_solid(material)
    k0 = 2π/λ0
    if prop
        grid = Grid.EnvGrid(thickness, λ0, (200e-9, 4e-6), 2000e-15)
        Et = sqrt.(gauss.(grid.t, fwhm=τfwhm))
        Eω = FFTW.fft(Et)
        Fields.prop_material!(Eω, grid, material, -thickness, λ0)
        P0int, _ = hquadrature(0, thickness) do ti
            Eωprop = Fields.prop_material(Eω, grid, material, ti, λ0)
            Etprop = FFTW.ifft(Eωprop)
            maximum(abs2.(Etprop))/maximum(abs2.(Et)) * peakpower
        end
    else
        P0int = peakpower * thickness
    end
    I0int = 2*P0int/(π*w0^2)
    n2 * k0 * I0int
end

end