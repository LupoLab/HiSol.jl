module Data
import Luna.Maths: derivative
import Luna.PhysData: wlfreq, sellmeier_gas

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