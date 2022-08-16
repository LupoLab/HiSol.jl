module Limits
import Luna.PhysData: density, pressure, n2_gas, γ3_gas, ε_0, c, ionisation_potential
import Luna.Ionisation: barrier_suppression
import Luna.Tools: field_to_intensity
import HISOL.Solitons: τfwhm_to_T0, T0P0, Δβwg, Δβρ
import HISOL.HCF: Aeff0, δ, fβ2, get_unm, ZDW

function critical_power(gas, pressure, λ)
    # Fibich and Gaeta, Optics Letters 25, 335 (2000)
    # Note here we use n_gas ~ 1, which is justified even for high pressure
    # (n_gas - 1 < 1e-3 typically)
    # The critical power is therefore strictly inversely proportional to density
    if pressure <= 0
        return Inf
    end
    n2 = n2_gas(gas, pressure)
    if n2 == 0
        return Inf
    end
    return 1.86225 * λ^2/(4π*n2)
end


function critical_density(gas, λ, τfwhm, energy)
    # Note here we use n_gas ~ 1, which is justified even for high pressure
    # (n_gas - 1 < 1e-3 typically)
    _, P0 = T0P0(τfwhm, energy)

    γ3 = γ3_gas(gas)
    n2 = 1.86225 * λ^2/(4*π*P0)
    return 4/3*ε_0*c/γ3*n2
end


critical_pressure(gas, args...) = pressure(gas, critical_density(gas, args...))


function critical_intensity(λ_target, gas, λ0; kwargs...)
    # Pcrit = Aeff * Icrit = Aeff0 * a² * Icrit
    # Pcrit0/ρ = Aeff0 * a² * Icrit
    # => Icrit = Pcrit0/Aeff0 / (a² * ρ)

    # Pcrit = 1.86225 * λ²/(4*π*n2)
    # n2 = 3/4*γ3/(ε0*c) * ρ
    # => Pcrit = 1.86225 * λ²*ε0*c/(3*π*γ3*ρ)
    # => Pcrit0 = 1.86225 * λ²*ε0*c/(3*π*γ3)

    ρasq = Δβwg(λ_target, λ0; kwargs...)/Δβρ(λ_target, gas, λ0)
    γ3 = γ3_gas(gas)
    Pcrit0 = 1.86225 * λ0^2*ε_0*c/(3*π*γ3)

    return Pcrit0/Aeff0(; kwargs...)/ρasq
end

function barrier_suppression_intensity(gas)
    Ip = ionisation_potential(gas)
    Esupp = barrier_suppression(Ip, 1)
    field_to_intensity(Esupp)
end

function Nmax_ion(λzd, gas, λ0, τfwhm; S_ion=10, kwargs...)
    # eq. (S16) in Supplementary, Travers et al., Nat. Phot. 13, 547 (2019)
    Isupp = barrier_suppression_intensity(gas)
    n2_0 = n2_gas(gas, 1)/density(gas, 1) # n2 at reference pressure
    T0 = τfwhm_to_T0(τfwhm)
    u_nm = get_unm(kwargs...)
    sqrt(T0^2*n2_0*Isupp*u_nm^2 / (S_ion*π*λ0*abs(δ(gas, λ0, λzd; kwargs...))*fβ2(gas, λzd)))
end

function Nmax_ion(a, gas, pressure, λ0, τfwhm; S_ion=10, kwargs...)
    Nmax_ion(ZDW(a, gas, pressure; kwargs...), gas, λ0, τfwhm; S_ion, kwargs...)
end

function Nmax_sf(λzd, gas, λ0, τfwhm; S_sf=5, kwargs...)
    # eq. (S15) in Supplementary, Travers et al., Nat. Phot. 13, 547 (2019)
    # But with factor in Pcrit of 1.86225 instead of 3
    T0 = τfwhm_to_T0(τfwhm)
    sqrt(1.86225T0^2*λ0/(3*S_sf*abs(δ(gas, λ0, λzd))))
end

function Nmax_sf(a, gas, pressure, λ0, τfwhm; S_sf=5, kwargs...)
    Nmax_sf(ZDW(a, gas, pressure; kwargs...), gas, λ0, τfwhm; S_sf, kwargs...)
end

function Nmax(args...; S_ion=10, S_sf=5, kwargs...)
    min(
        Nmax_ion(λzd, gas, λ0, τfwhm; S_ion, kwargs...),
        Nmax_sf(λzd, gas, λ0, τfwhm; S_sf, kwargs...)
    )
end

end