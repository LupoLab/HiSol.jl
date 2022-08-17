module Solitons
import HISOL: HCF, Data
import HISOL.HCF: intensity_modeavg
import Roots: find_zero
import Luna.PhysData: wlfreq, c, ε_0, γ3_gas, pressure

const T0fac = (2*acosh(sqrt(2)))

τfwhm_to_T0(τfwhm) = τfwhm / T0fac

T0_to_τfwhm(T0) = T0 * T0fac

function T0P0(τfwhm, energy)
    T0 = τfwhm_to_T0(τfwhm)
    P0 = energy/2T0
    T0, P0
end

function intensity_modeavg(a, τfwhm, energy; kwargs...)
    _, P0 = T0P0(τfwhm, energy)
    intensity_modeavg(a, P0; kwargs...)
end

function β_sol(a, gas, pressure, λ0, τfwhm, energy; kwargs...)
    _, P0 = T0P0(τfwhm, energy)
    γ_ = HCF.γ(a, gas, pressure, λ0; kwargs...)
    γ_*P0/2
end

function N(a, gas, pressure, λ0, τfwhm, energy; kwargs...)
    β2 = HCF.dispersion(a, gas, pressure, λ0, 2; kwargs...)
    T0, P0 = T0P0(τfwhm, energy)
    γ_ = HCF.γ(a, gas, pressure, λ0)
    sqrt(γ_*P0*T0^2/abs(β2))
end

function nonlinear_length(a, gas, pressure, λ0, τfwhm, energy; kwargs...)
    _, P0 = T0P0(τfwhm, energy)
    γ_ = HCF.γ(a, gas, pressure, λ0)
    1/(γ_*P0)
end

function dispersion_length(a, gas, pressure, λ0, τfwhm; kwargs...)
    T0 = τfwhm_to_T0(τfwhm)
    β2 = HCF.dispersion(a, gas, pressure, λ0, 2; kwargs...)
    T0^2/abs(β2)
end

function fission_length(a, gas, pressure, λ0, τfwhm, energy; kwargs...)
    L_d = dispersion_length(a, gas, pressure, λ0, τfwhm; kwargs...)
    L_nl = nonlinear_length(a, gas, pressure, λ0, τfwhm, energy; kwargs...)
    sqrt(L_d*L_nl)
end

function RDW_wavelength(a, gas, pressure, λ0, τfwhm=Inf, energy=0; kwargs...)
    ω0 = wlfreq(λ0)
    λzd = HCF.ZDW(a, gas, pressure; kwargs...)
    β1 = HCF.dispersion(a, gas, pressure, λ0, 1; kwargs...)
    β0 = HCF.β(a, gas, pressure, λ0; kwargs...)
    β_s = β_sol(a, gas, pressure, λ0, τfwhm, energy; kwargs...)

    lbω = wlfreq(λzd)
    ubω = wlfreq(Data.λmin[gas])

    ωRDW = missing
    try
        ωRDW = find_zero((lbω, ubω)) do ω
            HCF.β(a, gas, pressure, wlfreq(ω); kwargs...) - (ω-ω0)*β1 - β0 - β_s
        end
    catch e
        println(e)
    end
    return wlfreq(ωRDW)
end

function Δβwg(λ_target, λ0; kwargs...)
    u_nm = HCF.get_unm(kwargs...) # Bessel 0
    ωRDW = wlfreq(λ_target) # RDW frequency
    ω0 = wlfreq(λ0) # pump frequency
    u_nm^2*c/2*(1/ωRDW - 2/ω0 + ωRDW/ω0^2)
end

function Δβρ(λ_target, gas, λ0)
    ωRDW = wlfreq(λ_target) # RDW frequency
    ω0 = wlfreq(λ0) # pump frequency
    γ1_0 = Data.γ1(gas)(λ0) # γ⁽¹⁾ at pump freq
    γ1_RDW = Data.γ1(gas)(λ_target) # γ⁽¹⁾ at RDW freq
    γ11_0 = Data.dγ1dω(gas, λ0, 1) # dγ⁽¹⁾/dω at pump freq
    (ωRDW*γ1_RDW - ωRDW*γ1_0 - ωRDW*ω0*γ11_0 + ω0^2*γ11_0)/(2*c)
end

function RDW_density(λ_target, a, gas, λ0, τfwhm=Inf, energy=0; kwargs...)
    top = Δβwg(λ_target, λ0; kwargs...)/a^2
    bot = Δβρ(λ_target, gas, λ0)
    if energy > 0
        ~isfinite(τfwhm) && error("If energy is given, duration must also be given!")
        _, P0 = T0P0(τfwhm, energy)
        γ3 = γ3_gas(gas)
        aeff = HCF.Aeff(a; kwargs...)
        ω0 = wlfreq(λ0)
        bot -= 3*ω0*P0*γ3/(8*c^2*ε_0*aeff)
    end

    top/bot
end

RDW_pressure(λ_target, a, gas, args...; kwargs...) = pressure(
    gas,
    RDW_density(λ_target, a, gas, args...; kwargs...)
)

"""
    RDW_to_ZDW(λ0, λ_target, gas; kwargs...)

Given pump wavelength `λ0`, RDW wavelength `λ_target`, and the `gas` species,
find the corresponding zero-dispersion wavelength.
"""
function RDW_to_ZDW(λ0, λ_target, gas; kwargs...)
    ρasq_rdw = Δβwg(λ_target, λ0; kwargs...)/Δβρ(λ_target, gas, λ0)
    u_nm = HCF.get_unm(kwargs...)

    ω_target = wlfreq(λ_target)
    ω0 = wlfreq(λ0)
    ωguess = (ω_target + 2ω0)/3 # calculated value for only GVD (β₂) and TOD (β₃)

    fβ2 = HCF.fβ2(gas)
    ωzd = missing
    try
        ωzd = find_zero(ωguess) do ω
            u_nm^2/(2π^2*fβ2(wlfreq(ω))) - ρasq_rdw
        end
    catch
    end

    return wlfreq(ωzd)
end

end