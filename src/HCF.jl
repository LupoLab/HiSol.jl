module HCF
import Luna.PhysData: ref_index, c, ε_0, wlfreq, sellmeier_gas, density
import Luna.Capillary: get_unm, Aeff_Jintg
import Luna.Maths: derivative
import Roots: find_zero
import HISOL: Data

export α, β, dB_per_m, loss_length, transmission, Leff, dispersion, Aeff, Aeff0

get_unm(;n=1, m=1, kind=:HE) = get_unm(n, m, kind)

function α(a::Number, λ::Number; kwargs...)
    u_nm = get_unm(kwargs...)
    ω = wlfreq(λ)
    ν = real(ref_index(:SiO2, λ)) # ν = n_glass/n_gas with n_gas ≈ 1
    2*c^2*u_nm^2/(a^3*ω^2) * (ν^2 + 1)/(2*sqrt(ν^2-1))
end

function αbar_a(λ::Number; kwargs...)
    u_nm = get_unm(kwargs...)
    ω = wlfreq(λ)
    ν = real(ref_index(:SiO2, λ)) # ν = n_glass/n_gas with n_gas ≈ 1
    2*c^2*u_nm^2/ω^2 * (ν^2 + 1)/(2*sqrt(ν^2-1))
end

dB_per_m(args...; kw...) = 10/np.log(10)*α(args...; kw...)

loss_length(args...; kw...) = 1/α(args...; kw...)

transmission(length, a...; kw...) = exp(-α(a...; kw...)*length) 

function Leff(length, a...; kw...)
    α_ = α(a...; kw...)
    (1-exp(-α_*length))/α_
end

function βfunω(a, gas, pressure; kwargs...)
    u_nm = get_unm(kwargs...)
    γ1fun = Data.γ1(gas)
    ρ = density(gas, pressure)
    ω -> ω/c * (1 + ρ*γ1fun(wlfreq(ω))/2 - u_nm^2*c^2/(2a^2*ω^2))
end

function β(a, gas, pressure, λ::Number; kwargs...)
    u_nm = get_unm(kwargs...)
    ω = wlfreq(λ)
    γ1 = real(sellmeier_gas(gas)(1e6λ))
    ρ = density(gas, pressure)
    ω/c * (1 + ρ*γ1/2 - u_nm^2*c^2/(2a^2*ω^2))
end

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

Aeff(a; kwargs...) = a^2 * Aeff0(kwargs...)

intensity_modeavg(a, P0; kwargs...) = P0/Aeff(a; kwargs...)

function γ(a, gas, pressure, λ0; kwargs...)
    n2 = Data.n2_gas(gas, pressure)
    ω0 = wlfreq(λ0)
    ω0/c*n2/Aeff(a; kwargs...)
end

function γ(a, gas, λ0; λzd, kwargs...)
    n20 = Data.n2_0(gas)
    u_nm = get_unm(kwargs...)
    f = fβ2(gas, λzd)
    2*n20*u_nm^2/(3π*a^4*λ0*f)
end

# eq. S5 of Supplementary, Travers et al., Nat. Phot. 13, 547 (2019)
fβ2(gas, λ) = Data.dγ1dλ(gas, λ, 2)
fβ2(gas) = Data.dγ1dλ(gas, 2)

function δ(gas, λ, λzd; kwargs...)
    # eq. S7 of Supplementary, Travers et al., Nat. Phot. 13, 547 (2019)
    get_unm(kwargs...)^2 * λ^3/(8π^3*c^2) * (fβ2(gas, λ)/fβ2(gas, λzd) - 1)
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
    get_unm(kwargs...)^2/(2*π^2*a^2*fβ2(gas, λzd))
end

ZDW_pressure(λzd, a, gas; kwargs...) = pressure(gas, ZDW_density(λzd, a, gas; kwargs...))

end