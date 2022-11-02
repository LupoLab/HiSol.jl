module GasFlow
import HiSol.Data: mean_speed, gas_viscosity
import Roots: find_zero
import Luna.PhysData: atm, bar

"""
    pumps

Dictionary with some example pump speeds for small/medium/large turbos and
a roughing pump.
"""
pumps = Dict(
    "SmallTurbo" => 43e-3, # Pfeiffer HiPace 60P
    "MediumTurbo" => 655e-3, # Pfeiffer HiPace 600
    "LargeTurbo" => 2000e-3, # Pfeiffer HiPace 2300
    "RoughingPump" => 15/3600 # Edwards nXDS15i
)

function tube_conductance(a, flength, gas, P1, P2=0)
    # P1 and P2 have to be in SI units (Pascal)!
    ID = 2a
    c = mean_speed(gas)
    η = gas_viscosity(gas)
    mol = π/12 * c*ID^3/flength
    visc = π/256*1/η*ID^4/flength*(P1+P2)
    return mol+visc
end

function tube_PVflow(a, flength, gas, P1, P2=0)
    # P1 and P2 have to be in SI units (Pascal)!
    Pmax = max(P1, P2)
    Pmin = min(P1, P2)
    return tube_conductance(a, flength, gas, P1, P2)*(Pmax - Pmin)
end

function gradient_end_pressure(a, flength, gas, Pmax, pump_speed)
    # Pmax has to be SI units (Pascal)!
    find_zero(1e-3) do p2
        tube_PVflow(a, flength, gas, Pmax, p2) - p2*pump_speed
    end
end

Pascal_to_mbar(P_Pascal) = P_Pascal/100
slpm(qpv) = qpv/atm * 1000 * 60
mbarlps(qpv) = qpv / bar * 1000 * 1000

function choked_pressure(a, flength, gas, Pmax)
    ID = 2a
    η = gas_viscosity(gas)
    c = 1007 #speed of sound--TODO
    ID^2*Pmax^2/(64η*flength*c)
end

function tube_PVflow_choked(a, flength, gas, Pmax)
    c = 1007 # TODO
    pd = choked_pressure(a, flength, gas, Pmax)
    pd * π*a^2 * c
end

function gradient_end_pressure_choked(a, flength, gas, Pmax, pump_speed)
    qpv = tube_PVflow_choked(a, flength, gas, Pmax)
    qpv/pump_speed # assumes constant pump speed. not true!
end

end