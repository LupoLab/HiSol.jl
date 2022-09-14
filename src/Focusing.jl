module Focusing
import HiSol.Data: n2_solid
import HiSol.Solitons: T0P0

"""
    window_distance(a, λ0, energy, τfwhm, thickness; material=:SiO2, zr_frac=0.2)

Calculate the minimum distance between the exit of an HCF with radius `a` and a window made of `material`
with a given `thickness` such that the Kerr lens for a pulse at wavelength `λ0` with `energy` and duration `τfwhm`
does not cause a shift in the focal length by more than `zr_frac` times the Rayleigh length.

Note: A shift of `zr_frac` = 0.2 causes a drop in coupling efficiency of approximately 1%.
"""
function window_distance(a, λ0, energy, τfwhm, thickness; material=:SiO2, zr_frac=0.2)
    n2 = n2_solid(material)
    _, P0 = T0P0(τfwhm, energy)
    w0 = 0.64a
    zr = rayleigh(w0, λ0)

    sqrt(8*n2*thickness*P0/zr_frac/π/w0^4*zr^3) # approximate
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

- [`window_distance`](@ref): `thickness`, `material`, `zr_frac`
- [`mirror_distance`](@ref): `LIDT`, `S_fluence`
"""
function max_flength(a, λ0, energy, τfwhm, maxlength;
                     thickness=1e-3, material=:SiO2, zr_frac=0.2,
                     LIDT=2000, S_fluence=5,
                     entrance_window=true, exit_window=true)

    dwin = window_distance(a, λ0, energy, τfwhm, thickness; material, zr_frac)
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
diverged_beam(a, λ0, distance) = distance*λ0/(π*0.64a)

end