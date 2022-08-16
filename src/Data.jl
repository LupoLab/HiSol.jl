module Data
import Luna.Maths: derivative
import Luna.PhysData: wlfreq, sellmeier_gas

# minimum wavelengths to use in root-finding functions
# determined by first UV resonance in the sellmeier expansion
λmin = Dict(
    :Xe => 115e-9,
    :Kr => 102e-9,
    :Ar => 91e-9,
    :Ne => 77e-9,
    :He => 89e-9,
    :HeJ => 57e-9
)

function γ1(gas)
    s = sellmeier_gas(gas)
    λ -> s(1e6λ)
end

γ1(gas, λ) = γ1(gas)(λ)

function dγ1dω(gas, λ, n)
    ω = wlfreq(λ)
    f = γ1(gas) # γ1(λ)
    derivative(wlfreq(λ), n) do ω
        f(wlfreq(ω))
    end
end

dγ1dλ(gas, λ, n) = derivative(γ1(gas), λ, n)

end