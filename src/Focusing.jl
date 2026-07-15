module Focusing
import HiSol.Data: n2_solid, elastic_limit
import HiSol.Solitons: T0P0
import Luna.Maths: gauss, fwhm
import Luna: FFTW, Grid, Fields
import Luna.Modes: hquadrature
import Luna.PhysData: c
import LinearAlgebra: mul!, ldiv!
import Luna.Capillary: besselj, get_unm
import Roots: find_zero
import HiSol: Plotting
import Logging: NullLogger, with_logger

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

function window_thickness_nonlinear(a, λ0, energy, τfwhm, distance; shape=:sech, kwargs...)
    _, peakpower = T0P0(τfwhm, energy; shape)
    window_thickness_nonlinear(a, λ0, peakpower, distance; kwargs...)
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
- Eksma femtoline AR coating for 1030 nm: 0.1 J/cm² = 1000 J/m²
"""
function mirror_distance(a, λ0, energy, LIDT; S_fluence=5)
    w0 = 0.64a
    zr = rayleigh(w0, λ0)
    maxfluence = LIDT/S_fluence
    sqrt(2energy*zr^2/(π*w0^2*maxfluence))
end

"""
    max_flength(a, energy, τfwhm, maxlength, input_constraint, output_constraint; pressure=nothing)

Find the maximum HCF length which can fit into a space of `maxlength` for a core radius `a`
and a pulse with `energy` and FWHM duration `τfwhm`, given the [`LengthConstraint`](@ref)s
`input_constraint` and `output_constraint` on the entrance and exit side of the HCF.
`pressure` (in bar) is required by constraints which depend on it (e.g. [`WindowConstraint`](@ref)).
"""
function max_flength(a, energy, τfwhm, maxlength, input_constraint, output_constraint;
                     pressure=nothing)
    L_in = input_constraint(a, energy, τfwhm; pressure)
    L_out = output_constraint(a, energy, τfwhm; pressure)
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
        grid = with_logger(NullLogger()) do 
            Grid.EnvGrid(thickness, λ0, (100e-9, 4e-6), 2000e-15)
        end
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

function plot_window_thickness_variable(a, pressure, energy, τfwhm, λ0, λmax;
                            Bmax=0.2, n2=:SiO2, elastic_limit=:SiO2, S_break=4,
                            aperture_factor=2,
                            LIDT=nothing, S_fluence=5,
                            max_aperture_radius=10e-3)
    maxdist = beamsize_distance(0.64a, λmax, max_aperture_radius/aperture_factor)
    distance = collect(range(0, maxdist, 512))

    ΔP = max(pressure-1, 1) # make sure we can handle vacuum

    if ~isnothing(LIDT)
        LIDT_distance = mirror_distance(a, λ0, energy, LIDT; S_fluence)
    end

    w0win = diverged_beam.(a, λmax, distance)
    tNL = window_thickness_nonlinear.(a, λ0, energy, τfwhm, distance; material=n2, Bmax)
    tP = window_thickness_breaking.(ΔP, aperture_factor*w0win, elastic_limit; S_break)

    dOpt, tOpt, apOpt = window_distance_thickness_aperture(
        a, pressure, τfwhm, energy, λ0, λmax;
        Bmax, n2, elastic_limit, S_break, aperture_factor,
        max_aperture_radius=4*max_aperture_radius)

    plotdata = (;distance, tNL, tP,
                 aperture=aperture_factor*w0win,
                 dOpt, tOpt, apOpt,
                 LIDT_distance=isnothing(LIDT) ? nothing : LIDT_distance)
    Plotting.getext().plot_window_thickness_variable(plotdata)
end

function plot_window_thickness_fixed(a, pressure, peakpower, λ0, λmax, aperture_radius;
                                        Bmax=0.2, material=:SiO2, aperture_factor=2)
    mindist = rayleigh(0.64a, λmax)
    maxdist = beamsize_distance(0.64a, λmax, aperture_radius/aperture_factor)
    distance = collect(range(mindist, maxdist, 512))

    tNL = window_thickness_nonlinear.(a, λ0, peakpower, distance; material, Bmax)
    tP = window_thickness_breaking(pressure-1, aperture_radius, material)

    idx = findfirst(eachindex(distance)) do ii
        tNL[ii] >= tP
    end

    plotdata = (;distance, tNL, tP,
                 crossing=isnothing(idx) ? nothing : (;distance=distance[idx], thickness=tNL[idx]))
    Plotting.getext().plot_window_thickness_fixed(plotdata)
end

"""
    window_thickness_distance(a, pressure, τfwhm, energy, λ0, λmax, aperture_radius;
        Bmax=0.2, material=:SiO2, aperture_factor=2, round_thickness=true, S_break=4.0,
        shape=:sech)
    window_thickness_distance(a, pressure, peakpower, λ0, λmax, aperture_radius;
        Bmax=0.2, material=:SiO2, aperture_factor=2, round_thickness=true, S_break=4.0)

Calculate the required window thickness and the mimimum and maximum distance between the
exit of an HCF and the window.

Calculate the required window thickness to hold the specified `pressure` for the specified
window `material` and `aperture_radius`, and mimimum and maximum distance between the exit
of an HCF with radius `a` and the window, such that the B-integral for a pulse at
wavelength `λ0` with `energy` and duration `τfwhm` (or a given `peakpower`) does not exceed
`Bmax`, and such that the beam at wavelength `λmax` does not diverge such that the Gaussian
beam radius is less than `aperture_factor` smaller than `aperture_radius`.
`round_thickness` can be `true` (round to nearest mm, default), `false` (do not round)
or `x::Number` (round to the nearest x mm, e.g. `x=0.5` to round to nearest half mm)
The window thickness is calculated accounting for the safety factor `S_break`.
"""
function window_thickness_distance(a, pressure, τfwhm, energy, λ0, λmax, aperture_radius;
                                   shape=:sech, kwargs...)
    _, P0 = T0P0(τfwhm, energy; shape)
    window_thickness_distance(a, pressure, P0, λ0, λmax, aperture_radius; kwargs...)
end

function window_thickness_distance(a, pressure, peakpower, λ0, λmax, aperture_radius;
                                       Bmax=0.2, material=:SiO2, aperture_factor=2,
                                       round_thickness=true, S_break=4.0,
                                       elastic_limit=:SiO2)
    maximum_distance = beamsize_distance(0.64a, λmax, aperture_radius/aperture_factor)
    pressure = max(pressure - 1.0, 1.0) # make sure can handle vacuum
    thickness = window_thickness_breaking(pressure, aperture_radius, elastic_limit; S_break)
    if needround(round_thickness)
        rt = rounding(round_thickness)*1e-3
        thickness = ceil(thickness/rt)*rt
    end
    minimum_distance = window_distance(a, λ0, peakpower, thickness; material, Bmax)
    if minimum_distance > maximum_distance
        error("unable to find solution for given parameters")
    end
    (; thickness, minimum_distance, maximum_distance)
end

needround(r::Number) = true
needround(r::Bool) = r

rounding(r::Number) = r
rounding(r::Bool) = 1

"""
    window_distance_thickness_aperture(a, pressure, τfwhm, energy, λ0, λmax; kwargs...)

Find the window distance, thickness and aperture radius which satisfy both the nonlinearity constraint
and the pressure constraint with the shortest distance, thinnest window and smallest aperture possible.

# Arguments: 
- `a`: core radius
- `pressure`: gas pressure in bar
- `τfwhm`: pulse duration FWHM
- `energy`: pulse energy
- `λ0`: reference wavelength (for nonlinearity calculation)
- `λmax`: longest wavelength which needs to pass through the window unobstructed

# Keyword arguments:
- `Bmax`: maximum B-integral (default: 0.2)
- `n2`: nonlinear ref. index of the window, can be `Symbol` (material) or `Number` (n2)
    (default: :SiO2)
- `elastic_limit`: elastic limit of the window, can be `Symbol` (material) or `Number` (limit in Pascal)
    (default: :SiO2)
- `S_break`: safety factor for pressure tolerance calculation (default: 4)
- `aperture_factor`: multiple of 1/e² radius which determines the aperture radius (default: 2)
- `round_thickness`: whether to round the thickness to a nearest round value.
    Can be `true` (round to nearest mm), `false` (do not round)
    or `x::Number` (round to the nearest x mm, e.g. `x=0.5` to round to nearest half mm)
    (default: `false`)
- `round_aperture`: like `round_thickness`, but for aperture radius (default: `false`)
- `shape`: pulse shape to assume in calculating the peak power. can be :gauss or :sech
    (default: `sech`)
- `max_aperture_radius`: largest aperture radius to consider when searching for the optimal
    distance (default: 50 mm)
- `err`: if `true` (default), throw an error when no valid window exists below
    `max_aperture_radius`; if `false`, return `(Inf, Inf, Inf)` instead
"""
function window_distance_thickness_aperture(a, pressure, τfwhm, energy, λ0, λmax;
                                            Bmax=0.2, n2=:SiO2,
                                            elastic_limit=:SiO2, S_break=4,
                                            aperture_factor=2,
                                            round_thickness=false, round_aperture=false,
                                            shape=:sech,
                                            max_aperture_radius=50e-3,
                                            err=true)
    ΔP = max(pressure-1, 1) # make sure we can handle vacuum

    tdiff(d) = let w0win = diverged_beam(a, λmax, d)
        tNL = window_thickness_nonlinear(a, λ0, energy, τfwhm, d; material=n2, Bmax, shape)
        tP = window_thickness_breaking(ΔP, aperture_factor*w0win, elastic_limit; S_break)
        tNL - tP
    end

    if tdiff(0) > 0
        dOpt = 0
    else
        maxdist = beamsize_distance(0.64a, λmax, max_aperture_radius/aperture_factor)
        if tdiff(maxdist) < 0
            err && error(
                "Nonlinear and pressure thickness limits do not cross below an aperture "*
                "radius of $(1e3max_aperture_radius) mm; increase max_aperture_radius.")
            # no valid window exists below max_aperture_radius: signal an impossible
            # constraint (the convention of the other WindowConstraint branches) so that
            # design-space scans can treat this point as infeasible instead of crashing
            return Inf, Inf, Inf
        end
        dOpt = find_zero(tdiff, (0, maxdist))
    end

    distance = dOpt
    thickness = window_thickness_nonlinear(a, λ0, energy, τfwhm, dOpt;
                                           material=n2, shape, Bmax)
    aperture = aperture_factor*diverged_beam(a, λmax, dOpt)
    if needround(round_thickness) && needround(round_aperture)
        rt = rounding(round_thickness) * 1e-3
        ra = rounding(round_aperture) * 1e-3
        n = 0
        while true
            # round up aperture
            aperture = ceil(aperture/ra)*ra
            # need a thicker window for larger aperture
            tP = window_thickness_breaking(ΔP, aperture, elastic_limit; S_break)
            # round up thickness
            thickness = ceil(max(tP, thickness)/rt)*rt
            # distance required from nonlinearity constraint for new thickness
            distance = window_distance(a, λ0, energy, τfwhm, thickness;
                                       material=n2, shape, Bmax)
            # aperture required for new distance
            min_aperture = aperture_factor*diverged_beam(a, λmax, distance)
            n += 1
            if aperture >= min_aperture # all is well
                break
            elseif n >= 20
                error("Could not find window parameters")
            else
                aperture = min_aperture
            end
        end
    elseif needround(round_thickness)
        rt = rounding(round_thickness) * 1e-3
        n = 0
        while true
            # round up thickness
            thickness = ceil(thickness/rt)*rt
            # thicker window needs to be further away
            distance = window_distance(a, λ0, energy, τfwhm, thickness;
                                       material=n2, shape, Bmax)
            # need bigger aperture for larger distance
            aperture = aperture_factor*diverged_beam(a, λmax, distance)
            # thickness required for new aperture
            min_thickness = window_thickness_breaking(ΔP, aperture, elastic_limit; S_break)
            n += 1
            if thickness >= min_thickness # all is well
                break
            elseif n >= 20
                error("Could not find window parameters")
            else # new min thickness is too large -- need to go up one step further and try again
                thickness = min_thickness 
            end
        end
    elseif needround(round_aperture)
        ra = rounding(round_aperture) * 1e-3
        n = 0
        while true
            # round up aperture
            aperture = ceil(aperture/ra)*ra
            # may need thicker window for larger aperture
            thickness = window_thickness_breaking(ΔP, aperture, elastic_limit; S_break)
            # thicker window needs to be further away
            distance = window_distance(a, λ0, energy, τfwhm, thickness;
                                       material=n2, shape, Bmax)
            # aperture required for new distance
            min_aperture = aperture_factor*diverged_beam(a, λmax, distance)
            n += 1
            if aperture >= min_aperture # all is well
                break
            elseif n >= 20
                error("Could not find window parameters")
            else
                aperture = min_aperture
            end
        end
    end

    # Double check that our numbers are correct
    # 1. Is the window thin enough at this distance to avoid excessive nonlinearity?
    max_thickness = window_thickness_nonlinear(a, λ0, energy, τfwhm, distance;
                                               material=n2, shape, Bmax)
    if thickness > max_thickness + eps(Float64) # allow for floating-point differences
        error("Found thickness $(1e3thickness) mm which is too thick (max $(1e3max_thickness))")
    end
    # 2. Is the window thick enough to hold the pressure?
    min_thickness = window_thickness_breaking(ΔP, aperture, elastic_limit; S_break)
    if thickness < min_thickness - eps(Float64) # allow for floating-point differences
        error("Found thickness $(1e3thickness) mm which is too thin (min $(1e3min_thickness) mm )")
    end
    # 3. Is the aperture large enough to pass the beam through?
    min_aperture = aperture_factor*diverged_beam(a, λmax, distance)
    if aperture < min_aperture - eps(Float64) # allow for floating-point differences
        error("Found aperture $(1e3aperture) mm which is too small (min $(1e3min_aperture))")
    end
    return distance, thickness, aperture

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

"""
    LengthConstraint

Abstract supertype for length constraints on HCF systems. A `LengthConstraint` is *callable*
with the signature

    (c::LengthConstraint)(a, energy, τfwhm; pressure=nothing)

and returns the minimum distance between the HCF entrance/exit and whatever the constraint
refers to (window, mirror etc.) for core radius `a`, pulse `energy` and FWHM duration `τfwhm`.
`pressure` is the gas pressure in the HCF in bar; constraints which do not depend on it accept
it as an unused keyword argument.

All `LengthConstraint`s also implement [`details`](@ref), which returns more detailed
information about the geometry the constraint implies.
"""
abstract type LengthConstraint end

"""
    details(c::LengthConstraint, a, energy, τfwhm; pressure=nothing)

Return a `NamedTuple` with detailed information about the geometry implied by the constraint
`c` for core radius `a`, pulse `energy` and FWHM duration `τfwhm` at the given `pressure`.
The `NamedTuple` contains at least the field `distance`, which is identical to the result of
`c(a, energy, τfwhm; pressure)`.
"""
details(c::LengthConstraint, args...; kwargs...) = (;distance=c(args...; kwargs...))

"""
    NoConstraint()

A [`LengthConstraint`](@ref) which imposes no constraint at all, i.e. the HCF can take up
the whole available length. Calling it returns zero distance.
"""
struct NoConstraint <: LengthConstraint end

(nc::NoConstraint)(a, energy, τfwhm; pressure=nothing) = 0.0

details(nc::NoConstraint, a, energy, τfwhm; pressure=nothing) = (;distance=0.0)

"""
    FixedConstraint(distance)

A [`LengthConstraint`](@ref) which imposes a fixed `distance` between the HCF and the
first/last other optical element, regardless of all other parameters.
"""
struct FixedConstraint <: LengthConstraint
    distance::Float64
end

(fc::FixedConstraint)(a, energy, τfwhm; pressure=nothing) = fc.distance

details(fc::FixedConstraint, a, energy, τfwhm; pressure=nothing) = (;distance=fc.distance)

"""
    DamageConstraint(λref, LIDT; S_fluence=5, conversion=1)

A [`LengthConstraint`](@ref) which keeps the fluence at wavelength `λref` on the first/last
mirror below the damage threshold `LIDT` (**in SI units**, i.e. J/m²) divided by the safety
factor `S_fluence`. `conversion` is an energy conversion factor between the pulse energy in
the HCF and the pulse energy hitting the mirror (e.g. a frequency-conversion efficiency on
the exit side).

`details` for this constraint additionally returns `beam_radius`, the 1/e² beam radius at
`λref` on the mirror.

See also [`mirror_distance`](@ref).
"""
struct DamageConstraint <: LengthConstraint
    λref::Float64
    LIDT::Float64
    S_fluence::Float64
    conversion::Float64
end

function DamageConstraint(λref, LIDT; S_fluence=5, conversion=1)
    DamageConstraint(λref, LIDT, S_fluence, conversion)
end

function (dc::DamageConstraint)(a, energy, τfwhm; pressure=nothing)
    mirror_distance(a, dc.λref, energy*dc.conversion, dc.LIDT; S_fluence=dc.S_fluence)
end

function details(dc::DamageConstraint, a, energy, τfwhm; pressure=nothing)
    distance = dc(a, energy, τfwhm; pressure)
    beam_radius = diverged_beam(a, dc.λref, distance)
    (;distance, beam_radius)
end

"""
    WindowConstraint(λref, n2; kwargs...)

A [`LengthConstraint`](@ref) for a window at the entrance/exit of the HCF, made of a material
with nonlinear refractive index `n2` (either a `Symbol` for a known material, or a `Number`
in m²/W). The constraint keeps the B-integral for a pulse at wavelength `λref` below `Bmax`,
while requiring that the window is thick enough to hold the gas pressure and that the beam at
`λmax` passes through the aperture unobstructed. Thickness and aperture can each be fixed or
left variable (found automatically); if a fixed `thickness` cannot hold the pressure, or a
fixed `aperture` cannot pass the beam, calling the constraint returns `Inf`.

`details` for this constraint additionally returns `thickness` and `aperture`.

# Keyword arguments:
- `λmax`: longest wavelength which needs to pass through the window unobstructed
    (default: `λref`)
- `Bmax`: maximum B-integral in the window (default: 0.2)
- `thickness`: window thickness; `Number` (fixed) or `nothing` (variable, default)
- `round_thickness`: whether to round a variable thickness to a nearby round value. Can be
    `true` (round up to the next mm), `false` (do not round, default) or `x::Number`
    (round up to the next multiple of x mm, e.g. `x=0.5` for half-mm steps)
- `aperture`: aperture radius; `Number` (fixed) or `nothing` (variable, default)
- `aperture_factor`: multiple of the 1/e² beam radius at `λmax` which determines the required
    aperture radius (default: 2)
- `round_aperture`: like `round_thickness`, but for the aperture radius (default: `false`)
- `LIDT`: damage threshold of the window in J/m²; `Number` (taken into account) or `nothing`
    (ignore window damage, default)
- `S_fluence`: safety factor on `LIDT` (default: 5)
- `conversion`: energy conversion factor between the pulse energy in the HCF and the energy
    passing through the window (default: 1)
- `elastic_limit`: elastic limit of the window material; `Symbol` (material) or `Number`
    (limit in Pascal). Defaults to `n2` if that is a material `Symbol`, and to `:SiO2`
    otherwise.
- `S_break`: safety factor on the pressure handling (default: 4)
- `τfwhm`: overrides the pulse duration for the nonlinearity calculation if required
    (default: `nothing`, i.e. use the duration passed when calling the constraint)
"""
struct WindowConstraint{LT, mT, rtT, raT, tT, eT, τT, aT} <: LengthConstraint
    λref::Float64 # reference wavelength
    λmax::Float64 # maximum wavelength we want to pass through unobstructed
    n2::mT # Symbol (material) or Number (n₂)
    Bmax::Float64 # Maximum B-integral
    thickness::tT # Number (fixed) or nothing (variable)
    round_thickness::rtT # true (round to next mm), false (do not round), or number (fraction of mm)
    aperture::aT # Fixed aperture (Number) or variable aperture (nothing)
    aperture_factor::Float64 # ratio between aperture radius and w₀ (1/e² radius)
    round_aperture::raT # true (round to next mm), false (do not round), or number (fraction of mm)
    LIDT::LT # Number (take into account window damage) or nothing (ignore window damage)
    S_fluence::Float64 # safety factor on LIDT
    conversion::Float64 # conversion factor between input and output energy
    elastic_limit::eT # Symbol (material) or Number (modulus of rupture)
    S_break::Float64 # Safety factor on pressure handling
    τfwhm::τT # Can overwrite pulse duration for nonlinearity calculation if required
end

function WindowConstraint(λref, n2;
                          λmax=λref, Bmax=0.2, thickness=nothing, round_thickness=false,
                          aperture=nothing, aperture_factor=2,
                          round_aperture=false,
                          LIDT=nothing, S_fluence=5,
                          conversion=1,
                          elastic_limit=nothing, S_break=4,
                          τfwhm=nothing)
    WindowConstraint(λref, λmax, n2, Bmax, thickness, round_thickness,
                     aperture, float(aperture_factor), round_aperture,
                     LIDT, float(S_fluence),
                     float(conversion), def_el(elastic_limit, n2), float(S_break),
                     τfwhm)
end

# find default elastic limit value
def_el(el::Nothing, n2::Symbol) = n2 # if n2 is given as a material, use that
def_el(el::Nothing, n2::Number) = :SiO2 # if n2 is given as number, use silica default
def_el(el, n2) = el # if elastic limit is given directly, use that

function (wc::WindowConstraint)(a, energy, τfwhm; pressure)
    details(wc, a, energy, τfwhm; pressure).distance
end

# Maximum distance between HCF and window such that the beam at λmax still fits through
# a fixed aperture; -Inf if it does not fit even at zero distance.
function _max_aperture_distance(a, λmax, aperture, aperture_factor)
    w0 = 0.64a
    if aperture > aperture_factor*w0
        beamsize_distance(w0, λmax, aperture/aperture_factor)
    else
        -Inf
    end
end

function details(wc::WindowConstraint, a, energy, τfwhm; pressure)
    energy *= wc.conversion
    τfwhm = isnothing(wc.τfwhm) ? τfwhm : wc.τfwhm
    ΔP = max(pressure-1, 1) # make sure we can handle vacuum
    if isnothing(wc.aperture) && isnothing(wc.thickness)
        # thickness and aperture are variable: find the joint optimum
        distance, thickness, aperture = window_distance_thickness_aperture(
            a, pressure, τfwhm, energy, wc.λref, wc.λmax;
            Bmax=wc.Bmax, n2=wc.n2,
            elastic_limit=wc.elastic_limit, S_break=wc.S_break,
            aperture_factor=wc.aperture_factor,
            round_thickness=wc.round_thickness,
            round_aperture=wc.round_aperture,
            err=false,
        )
    elseif isnothing(wc.aperture)
        # aperture is variable, thickness is fixed
        thickness = wc.thickness
        distance = window_distance(a, wc.λref, energy, τfwhm, thickness;
                                   material=wc.n2, Bmax=wc.Bmax)
        aperture = wc.aperture_factor * diverged_beam(a, wc.λmax, distance)
        if needround(wc.round_aperture)
            ra = rounding(wc.round_aperture)*1e-3
            aperture = ceil(aperture/ra)*ra
        end
        min_thickness = window_thickness_breaking(ΔP, aperture, wc.elastic_limit;
                                                  S_break=wc.S_break)
        if thickness < min_thickness
            # fixed thickness cannot hold the pressure at the required aperture
            distance = Inf
        end
    elseif isnothing(wc.thickness)
        # thickness is variable, aperture is fixed
        aperture = wc.aperture
        maximum_distance = _max_aperture_distance(a, wc.λmax, aperture, wc.aperture_factor)
        thickness = window_thickness_breaking(ΔP, aperture, wc.elastic_limit;
                                              S_break=wc.S_break)
        if needround(wc.round_thickness)
            rt = rounding(wc.round_thickness)*1e-3
            thickness = ceil(thickness/rt)*rt
        end
        distance = window_distance(a, wc.λref, energy, τfwhm, thickness;
                                   material=wc.n2, Bmax=wc.Bmax)
        if distance > maximum_distance
            # beam at λmax no longer fits through the fixed aperture
            distance = Inf
        end
    else
        # both thickness and aperture are fixed
        thickness = wc.thickness
        aperture = wc.aperture
        maximum_distance = _max_aperture_distance(a, wc.λmax, aperture, wc.aperture_factor)
        min_thickness = window_thickness_breaking(ΔP, aperture, wc.elastic_limit;
                                                  S_break=wc.S_break)
        distance = window_distance(a, wc.λref, energy, τfwhm, thickness;
                                   material=wc.n2, Bmax=wc.Bmax)
        if (thickness < min_thickness) || (distance > maximum_distance)
            distance = Inf
        end
    end
    if ~isnothing(wc.LIDT)
        dLIDT = mirror_distance(a, wc.λref, energy, wc.LIDT; S_fluence=wc.S_fluence)
        distance = max(distance, dLIDT)
    end
    (;distance, thickness, aperture)
end

end