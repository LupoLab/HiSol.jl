module Solitons
import HighEnergySolitons: CapillaryApproximations

T0fac = (2*acosh(sqrt(2)))

τfwhm_to_T0(τfwhm) = τfwhm / T0fac

T0_to_τfwhm(T0) = T0 * T0fac

function T0P0(τfwhm, energy)
    T0 = τfwhm_to_T0(τfwhm)
    P0 = energy/2T0
    T0, P0
end

function N(a, gas, pressure, λ0, τfwhm, energy; kwargs...)
    β2 = CapillaryApproximations.dispersion(a, λ0, gas, pressure, 2; kwargs...)
    T0, P0 = T0P0(τfwhm, energy)
    γ_ = CapillaryApproximations.γ(a, gas, pressure, λ0)
    sqrt(γ_*P0*T0^2/abs(β2))
end

function nonlinear_length(a, gas, pressure, λ0, τfwhm, energy; kwargs...)
    _, P0 = T0P0(τfwhm, energy)
    γ_ = CapillaryApproximations.γ(a, gas, pressure, λ0)
    1/(γ_*P0)
end

function dispersion_length(a, gas, pressure, λ0, τfwhm; kwargs...)
    T0 = τfwhm_to_T0(τfwhm)
    β2 = CapillaryApproximations.dispersion(a, λ0, gas, pressure, 2; kwargs...)
    T0^2/abs(β2)
end

function fission_length(a, gas, pressure, λ0, τfwhm, energy; kwargs...)
    L_d = dispersion_length(a, gas, pressure, λ0, τfwhm; kwargs...)
    L_nl = nonlinear_length(a, gas, pressure, λ0, τfwhm, energy; kwargs...)
    sqrt(L_d*L_nl)
end

end