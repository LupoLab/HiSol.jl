module HCF
import Luna.PhysData: ref_index, c, ε_0, wlfreq, sellmeier_gas, density, n2_gas
import Luna.Capillary: get_unm, Aeff_Jintg
import Luna.Maths: derivative
import Roots: find_zero

export α, β, dB_per_m, loss_length, transmission, Leff, dispersion, Aeff, Aeff0

get_unm(;n=1, m=1, kind=:HE) = get_unm(n, m, kind)

function α(a::Number, λ::Number; kwargs...)
    u_nm = get_unm(kwargs...)
    ω = wlfreq(λ)
    ν = ref_index(:SiO2, λ) # ν = n_glass/n_gas with n_gas ≈ 1
    2*c^2*u_nm^2/(a^3*ω^2) * (ν^2 + 1)/(2*sqrt(ν^2-1))
end

dB_per_m(a...; kw...) = 10/np.log(10)*α(a...; kw...)

loss_length(a...; kw...) = 1/α(args...; kw...)

transmission(length, a...; kw...) = exp(-α(a...; kw...)*length) 

function Leff(length, a...; kw...)
    α_ = α(a...; kw...)
    (1-exp(-α_*length))/α_
end

function β(a, λ::Number, gas, pressure; kwargs...)
    u_nm = get_unm(kwargs...)
    ω = wlfreq(λ)
    γ1 = real(sellmeier_gas(gas)(1e6λ))
    ρ = density(gas, pressure)
    ω/c * (1 + ρ*γ1/2 - u_nm^2*c^2/(2a^2*ω^2))
end

function dispersion(a, λ, gas, pressure, n=2; kwargs...)
    derivative(wlfreq(λ), n) do ω
        β(a, wlfreq(ω), gas, pressure; kwargs...)
    end
end

Aeff0(;n=1, m=1, kind=:HE) = Aeff_Jintg(n, get_unm(n, m, kind), kind)

Aeff(a; kwargs...) = a^2 * Aeff0(kwargs...)

function γ(a, gas, pressure, λ0; kwargs...)
    n2 = n2_gas(gas, pressure)
    ω0 = wlfreq(λ0)
    ω0/c*n2/Aeff(a; kwargs...)
end

# eq. S5 of Supplementary, Travers et al., Nat. Phot. 13, 547 (2019)
fβ2(gas, λ) = Data.dγ1dλ(gas, λ, 2)

function δ(gas, λ, λzd; kwargs...)
    # eq. S7 of Supplementary, Travers et al., Nat. Phot. 13, 547 (2019)
    get_unm(kwargs...)^2 * λ^3/(8π^3*c^2) * (fβ2(gas, λ)/fβ2(gas, λzd) - 1)
end

function ZDW(a, gas, pressure; λmin=100e-9, λmax=3e-6, kwargs...)
    ubω = 2π*c/λmin
    lbω = 2π*c/λmax
    ω0 = missing
    try
        ω0 = find_zero((lbω, ubω)) do ω
            dispersion(a, wlfreq(ω), gas, pressure, 2; kwargs...)
        end
    catch e
        println(e)
    end
    return 2π*c/ω0
end

end