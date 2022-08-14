module Limits
import Luna.PhysData: pressure, n2_gas, γ3_gas, ε_0, c
import HISOL.Solitons: T0P0, Δβwg, Δβρ
import HISOL.HCF: Aeff0

function critical_power(gas, pressure, λ)
    # Note here we use n_gas ~ 1, which is justified even for high pressure
    # (n_gas - 1 < 1e-3 typically)
    # The critical power is therefore strictly inversely proportional to density
    if pressure <= 0:
        return Inf
    end
    n2 = n2_gas(gas, pressure)
    if n2 == 0
        return Inf
    end
    return 1.86225 * λ**2/(4π*n2)
end


function critical_density(gas, λ, τfwhm, energy)
    # Note here we use n_gas ~ 1, which is justified even for high pressure
    # (n_gas - 1 < 1e-3 typically)
    _, P0 = T0P0(τfwhm, energy)

    γ3 = γ3_gas(gas)
    n2 = 1.86225 * λ**2/(4*π*P0)
    return 4/3*ε_0*c/γ3*n2
end


critical_pressure(gas, args...) = pressure(gas, critical_density(gas, args...))


function critical_intensity(λ_target, gas, λ0; kwargs...):
    # Pcrit = Aeff * Icrit = Aeff0 * a**2 * Icrit
    # Pcrit0/ρ = Aeff0 * a**2 * Icrit
    # => Icrit = Pcrit0/Aeff0 / (a**2 * ρ)

    # Pcrit = 1.86225 * λ**2/(4*π*n2)
    # n2 = 3/4*γ3/(ε0*c) * ρ
    # => Pcrit = 1.86225 * λ**2*ε0*c/(3*π*γ3*ρ)
    # => Pcrit0 = 1.86225 * λ**2*ε0*c/(3*π*γ3)

    ρasq = Δβwg(λ_target, λ0; kwargs...)/Δβρ(λ_target, gas, λ0)
    γ3 = γ3_gas(gas)
    Pcrit0 = 1.86225 * λ0^2*ε_0*c/(3*π*γ3)

    return Pcrit0/Aeff0(; kwargs...)/ρasq
end

end