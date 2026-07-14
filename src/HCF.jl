module HCF
import Luna.PhysData: ref_index, c, ε_0, wlfreq, sellmeier_gas, density, pressure
import Luna.Capillary: get_unm, Aeff_Jintg, MarcatiliMode
import Luna.Maths: derivative
import Luna.Fields: normalised_gauss_beam
import Luna.Modes: overlap
import Roots: find_zero
import HiSol: Data

export α, β, dB_per_m, loss_length, transmission, Leff, dispersion, Aeff, Aeff0

get_unm(;n=1, m=1, kind=:HE) = get_unm(n, m, kind)

"""
    α(a::Number, λ::Number; n=1, m=1, kind=:HE)

Calculate the power loss coefficient for an HCF with core radius `a` at
wavelength `λ` for the mode defined by `n`, `m`, and `kind`.
"""
α(a, args...; kwargs...) = αbar_a(args...; kwargs...)/a^3

"""
    αbar_a(λ::Number; n=1, m=1, kind=:HE)

Calculate the core-size independent prefactor of the power loss coefficient
for the HCF mode defined by `n`, `m`, and `kind`at wavelength `λ`.
"""
function αbar_a(λ::Number; n=1, m=1, kind=:HE)
    u_nm = get_unm(;n, m, kind)
    ω = wlfreq(λ)
    ν = real(ref_index(:SiO2, λ)) # ν = n_glass/n_gas with n_gas ≈ 1
    if kind == :HE
        return 2*c^2*u_nm^2/ω^2 * (ν^2 + 1)/(2*sqrt(ν^2-1))
    elseif kind == :TE
        return 2*c^2*u_nm^2/ω^2 * 1/sqrt(ν^2-1)
    elseif kind == :TM
        return 2*c^2*u_nm^2/ω^2 * ν^2/sqrt(ν^2-1)
    else
        error("Unknown mode kind $kind")
    end
end

"""
    dB_per_m(a::Number, λ::Number; n=1, m=1, kind=:HE)

Calculate the power loss coefficient in dB/m for an HCF with
core radius `a` at wavelength `λ` for the mode defined by
`n`, `m`, and `kind`.
"""
dB_per_m(args...; kw...) = 10/log(10)*α(args...; kw...)

"""
    loss_length(a::Number, λ::Number; n=1, m=1, kind=:HE)

Calculate the 1/e power loss length for an HCF with core radius `a` at
wavelength `λ` for the mode defined by `n`, `m`, and `kind`.
"""
loss_length(args...; kw...) = 1/α(args...; kw...)

transmission(length, a...; kw...) = exp(-α(a...; kw...)*length)

function transmission(length, a, λ, Nmodes=6, afrac=0.64)
    k = 2π/λ
    w0 = afrac*a
    Erθ = normalised_gauss_beam(k, w0)
    t = 0.0
    for m in 1:Nmodes
        mode = MarcatiliMode(a; m)
        η = abs2.(overlap(mode, Erθ))
        α_ = α(a, λ; m)
        t += η * exp(-α_*length)
    end
    t
end

function Leff(length, a...; kw...)
    α_ = α(a...; kw...)
    (1-exp(-α_*length))/α_
end

function βfunω(a, gas, pressure; kwargs...)
    u_nm = get_unm(;kwargs...)
    γ1fun = Data.γ1(gas)
    ρ = density(gas, pressure)
    ω -> ω/c * (1 + ρ*γ1fun(wlfreq(ω))/2 - u_nm^2*c^2/(2a^2*ω^2))
end

β(a, gas, pressure, λ::Number; kwargs...) = βfunω(a, gas, pressure; kwargs...)(wlfreq(λ))

function β_ret(a, gas, pressure, λ, λ0; kwargs...)
    ω = wlfreq(λ)
    ω0 = wlfreq(λ0)
    β0 = β(a, gas, pressure, λ0; kwargs...)
    β1 = dispersion(a, gas, pressure, λ0, 1; kwargs...)
    βλ = β.(a, gas, pressure, λ; kwargs...)
    @. βλ - β1*(ω - ω0) - β0
end

function dispersion(a, gas, pressure, λ, n=2; kwargs...)
    βf = βfunω(a, gas, pressure; kwargs...)
    derivative(βf, wlfreq(λ), n)
end

Aeff0(;n=1, m=1, kind=:HE) = Aeff_Jintg(n, get_unm(n, m, kind), kind)

Aeff(a; kwargs...) = a^2 * Aeff0(;kwargs...)

intensity_modeavg(a, P0; kwargs...) = P0/Aeff(a; kwargs...)

function γ(a, gas, pressure, λ0; kwargs...)
    n2 = Data.n2_gas(gas, pressure)
    ω0 = wlfreq(λ0)
    ω0/c*n2/Aeff(a; kwargs...)
end

function γ(gas, λ0; λzd=nothing, ρasq=nothing, kwargs...)
    n20 = Data.n2_0(gas)
    if ~isnothing(λzd)
        u_nm = get_unm(;kwargs...)
        f = fβ2(gas, λzd)
        return n20*u_nm^2/(Aeff0(;kwargs...) * π*λ0*f)
    elseif ~isnothing(ρasq)
        isnothing(λzd) || error("Only one of ρasq or λzd kwargs can be given")
        return 2π/λ0 * n20 * ρasq/Aeff0(;kwargs...)
    else
        error("One of ρasq or λzd kwargs must be given")
    end
end

γ(a, args...; kwargs...) = γ(args...; kwargs...) / a^4

# eq. S5 of Supplementary, Travers et al., Nat. Phot. 13, 547 (2019)
fβ2(gas, λ) = Data.dγ1dλ(gas, λ, 2)
fβ2(gas) = Data.dγ1dλ(gas, 2)

function δ(gas, λ, λzd; kwargs...)
    # eq. S7 of Supplementary, Travers et al., Nat. Phot. 13, 547 (2019)
    get_unm(;kwargs...)^2 * λ^3/(8π^3*c^2) * (fβ2(gas, λ)/fβ2(gas, λzd) - 1)
end

function Δ(gas, λ, ρasq; kwargs...)
    u_nm = get_unm(;kwargs...)
    ω = wlfreq(λ)
    ρasq/(2c) * (2*Data.dγ1dω(gas, λ, 1) + ω*Data.dγ1dω(gas, λ, 2)) - u_nm^2*c/ω^3
end

function ZDW(a, gas, pressure; λmin=Data.λmin[gas], λmax=3e-6, kwargs...)
    ubω = 2π*c/λmin
    lbω = 2π*c/λmax
    ω0 = missing
    try
        ω0 = find_zero((lbω, ubω)) do ω
            dispersion(a, gas, pressure, wlfreq(ω), 2; kwargs...)
        end
    catch e
        println(e)
    end
    return 2π*c/ω0
end

function ZDW_density(λzd, a, gas; kwargs...)
    # eq. S4 of Supplementary, Travers et al., Nat. Phot. 13, 547 (2019)
    get_unm(;kwargs...)^2/(2*π^2*a^2*fβ2(gas, λzd))
end

function ZDW_density_area_product(λzd, gas; kwargs...)
    get_unm(;kwargs...)^2/(2*π^2*fβ2(gas, λzd))
end

ZDW_pressure(λzd, a, gas; kwargs...) = pressure(gas, ZDW_density(λzd, a, gas; kwargs...))

end
