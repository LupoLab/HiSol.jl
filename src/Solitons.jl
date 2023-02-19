module Solitons
import HiSol: HCF, Data
import HiSol.HCF: intensity_modeavg
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

nonlinear_length(P0, γ_) = 1/(γ_*P0) 

function dispersion_length(a, gas, pressure, λ0, τfwhm; kwargs...)
    T0 = τfwhm_to_T0(τfwhm)
    β2 = HCF.dispersion(a, gas, pressure, λ0, 2; kwargs...)
    T0^2/abs(β2)
end

dispersion_length(T0, β2) = T0^2/abs(β2)

"""
    fission_length(a, gas, pressure, λ0, τfwhm, energy; kwargs...)

Calculate the fission length for HCF with core radius `a` filled with
`pressure` bars of `gas` when pumping at `λ0` with duration `τfwhm`
and pulse `energy`. Further `kwargs` `m`, `n` and `kind` determine the HCF mode.
"""
function fission_length(a, gas, pressure, λ0, τfwhm, energy; kwargs...)
    L_d = dispersion_length(a, gas, pressure, λ0, τfwhm; kwargs...)
    L_nl = nonlinear_length(a, gas, pressure, λ0, τfwhm, energy; kwargs...)
    sqrt(L_d*L_nl)
end

"""
    fission_length(a, gas, λ0, τfwhm; N, λzd, kwargs...)

Calculate the fission length for HCF with core radius `a` filled with `gas`
such that the ZDW is `λzd` when pumping at `λ0` with a duration `τfwhm` and
soliton order `N`. Further `kwargs` `m`, `n` and `kind` determine the HCF mode.
"""
function fission_length(a, gas, λ0, τfwhm; N, λzd, kwargs...)
    T0 = τfwhm_to_T0(τfwhm)
    δ = HCF.δ(gas, λ0, λzd; kwargs...)
    T0^2*a^2/(N*abs(δ))
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

function density_area_product(λ_target, gas, λ0; kwargs...)
    Δβwg(λ_target, λ0; kwargs...)/Δβρ(λ_target, gas, λ0; kwargs...)
end

"""
    RDW_density(λ_target, a, gas, λ0, τfwhm=Inf, energy=0; kwargs...)

Calculate the phase-matching `gas` density for RDW emission at `λ_target` in an
HCF with core radius `a` when pumping at `λ0`. The nonlinear contribution is taken
into account when both `τfwhm` and `energy` are given. The HCF mode can be
set by the keyword arguments `n`, `m`, and `kind`.
"""
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

"""
    RDW_pressure(λ_target, a, gas, λ0, τfwhm=Inf, energy=0; kwargs...)

Calculate the phase-matching `gas` pressure for RDW emission at `λ_target` in an
HCF with core radius `a` when pumping at `λ0`. The nonlinear contribution is taken
into account when both `τfwhm` and `energy` are given. The HCF mode can be
set by the keyword arguments `n`, `m`, and `kind`.
"""
RDW_pressure(λ_target, a, gas, args...; kwargs...) = pressure(
    gas,
    RDW_density(λ_target, a, gas, args...; kwargs...)
)

"""
    RDW_to_ZDW(λ0, λ_target, gas; kwargs...)

Given pump wavelength `λ0`, RDW wavelength `λ_target`, and the `gas` species,
find the corresponding zero-dispersion wavelength. Further `kwargs`
`m`, `n` and `kind` determine the HCF mode.
"""
function RDW_to_ZDW(λ0, λ_target, gas; kwargs...)
    ρasq_rdw = Δβwg(λ_target, λ0; kwargs...)/Δβρ(λ_target, gas, λ0)
    u_nm = HCF.get_unm(kwargs...)

    ω_target = wlfreq(λ_target)
    ω0 = wlfreq(λ0)
    ωguess = (ω_target + 2ω0)/3 # calculated value for only GVD (β₂) and TOD (β₃)

    fβ2 = HCF.fβ2(gas)
    ωzd = missing
    # fast method: start with initial guess
    try
        ωzd = find_zero(ωguess) do ω
            u_nm^2/(2π^2*fβ2(wlfreq(ω))) - ρasq_rdw
        end
    catch
    end

    # if fast method fails: λzd must be between λ_target and λ0
    # this is much slower but more reliable
    if ismissing(ωzd) || ωzd == 0
        try
            ωzd = find_zero((wlfreq(λ0), wlfreq(λ_target))) do ω
                u_nm^2/(2π^2*fβ2(wlfreq(ω))) - ρasq_rdw
            end
        catch
        end
    end

    return wlfreq(ωzd)
end

"""
    N_to_energy(N, a, gas, λ0, λzd, τfwhm; kwargs...)
    N_to_energy(N, a, gas, λ0, τfwhm; ρasq, kwargs...)

Convert soliton order `N` to energy for core radius `a` filled with `gas`
such that the ZDW is `λzd` and when pumping at `λ0` with a duration `τfwhm`.
Further `kwargs` `m`, `n` and `kind` determine the HCF mode. `ρasq` can be
given as a keyword argument instead of `λzd`.
"""
N_to_energy(N, a, args...; kwargs...) = N_to_energy(args...; kwargs...)(N, a)

function N_to_energy(gas::Symbol, λ0, λzd, τfwhm; kwargs...)
    T0 = τfwhm_to_T0(τfwhm)
    δ = HCF.δ(gas, λ0, λzd; kwargs...)
    γ = HCF.γ(gas, λ0; λzd, kwargs...)
    (N, a) -> 2*T0*N^2*abs(δ)/(T0^2*γ)*a^2
end

function N_to_energy(gas::Symbol, λ0, τfwhm; ρasq, kwargs...)
    T0 = τfwhm_to_T0(τfwhm)
    Δ = HCF.Δ(gas, λ0, ρasq; kwargs...)
    γ = HCF.γ(gas, λ0; ρasq, kwargs...)
    (N, a) -> 2*T0*N^2*abs(Δ)/(T0^2*γ)*a^2
end

end