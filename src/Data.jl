module Data
import Luna.Maths: derivative
import Luna.PhysData: wlfreq, sellmeier_gas, m_u, roomtemp, k_B, c, ε_0, γ3_gas, density, pressure

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

function dγ1dω(gas, n)
    f = γ1(gas) # γ1(λ)
    λ -> derivative(wlfreq(λ), n) do ω
        f(wlfreq(ω))
    end
end

dγ1dω(gas, λ, n) = dγ1dω(gas, n)(λ)

function dγ1dλ(gas, n)
    f = γ1(gas)
    λ -> derivative(f, λ, n)
end

dγ1dλ(gas, λ, n) = dγ1dλ(gas, n)(λ)

n2_0(gas) = 3/4 * γ3_gas(gas) / (ε_0*c)
n2_gas(gas, pressure) = n2_0(gas) * density(gas, pressure)


function n2_solid(material)
    if material == :SiO2
        return 2.6e-20
    elseif material == :MgF2
        # R. DeSalvo et al., IEEE J. Q. Elec. 32, 10 (1996).
        return 5.79e-21
    elseif material == :CaF2
        # Ensley Optics Express Vol. 27, No. 26 / 23 December 2019 
        return 1.8e-20
    elseif material == :CaCO3
        # Kabaciński et al., 10.1364/OE.27.011018
        return 3.22e-20
    else
        error("Unkown glass $material")
    end
end

function elastic_limit(material)
    # https://www.crystran.co.uk/media/wysiwyg/Datasheets/design-of-pressure-windows_1.pdf
    if material == :SiO2
        return 55e6
    elseif material == :MgF2
        return 49.6e6
    elseif material == :CaF2
        return 36.5e6
    elseif material == :Al2O3
        return 276e6
    else
        error("Unkown glass $material")
    end
end

function gas_viscosity(gas)
    if gas in (:He, :HeJ)
        return 19.6e-6
    elseif gas == :Ne
        return 31.9e-6
    elseif gas == :Ar
        return 22.11e-6
    elseif gas == :Kr
        return 1.4257*gas_viscosity(N2)
    elseif gas == :Xe
        return 1.2964*gas_viscosity(N2)
    elseif gas == :N2
        return 17.5e-6
    else
        error("Unknown gas $gas")
    end
end

function gas_mass(gas)
    if gas in (:He, :HeJ)
        m = 4
    elseif gas == :Ne
        m = 20.1797
    elseif gas == :Ar
        m = 39.948
    elseif gas == :Kr
        m = 83.798
    elseif gas == :Xe
        m = 131.293
    elseif gas == :H2
        m = 2
    elseif gas == :N2
        m = 14.0067
    else
        error("Unknown gas $gas")
    end

    return m*m_u
end

mean_speed(gas, T=roomtemp) = sqrt(8*k_B*T/(π*gas_mass(gas)))

function adiabatic_ratio(gas)
    if gas in (:He, :HeJ)
        return 1.66
    elseif gas == :Ne
        return 1.64
    elseif gas == :Ar
        return 1.67
    elseif gas == :Kr
        return 1.68
    elseif gas == :Xe
        return 1.66
    elseif gas == :H2
        return 1.41
    elseif gas == :N2
        return 1.47
    else
        error("Unknown gas $gas")
    end
end

function speed_of_sound(gas, T=roomtemp)
    γ = adiabatic_ratio(gas)
    m = gas_mass(gas)
    sqrt(γ*k_B*T/m)
end
end