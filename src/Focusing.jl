module Focusing
import HISOL.Data: n2_solid
import HISOL.Solitons: T0P0

"""
    window_distance(a, λ0, energy, τfwhm, thickness; material=:SiO2, zr_frac=0.2)

Calculate the minimum distance between the exit of an HCF with radius `a` and a window made of `material`
with a given `thickness` such that the Kerr lens for a pulse at wavelength `λ0` with `energy` and duration `τfwhm`
does not cause a shift in the focal length by more than `zr_frac` times the Rayleigh length.

Note: A shift of `zr_frac` = 0.2 causes a drop in coupling efficiency of approximately 1%.
"""
function window_distance(a, λ0, energy, τfwhm, thickness; material=:SiO2, zr_frac=0.5)
    n2 = n2_solid(material)
    _, P0 = T0P0(τfwhm, energy)
    w0 = 0.64a
    zr = rayleigh(w0, λ0)

    sqrt(8*n2*thickness*P0/zr_frac/π/w0^4*zr^3) # approximate
end

rayleigh(w0, λ) = π*w0^2/λ

end