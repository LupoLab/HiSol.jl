module GasFlow
import HiSol.Data: mean_speed, gas_viscosity, speed_of_sound
import Roots: find_zero
import Luna.PhysData: atm, bar
import Luna.Maths: BSpline
import CSV

datafolder = joinpath(dirname(dirname(@__FILE__)), "data")

function get_characteristic(filename)
    f = CSV.File(joinpath(datafolder, filename))
    pcol, scol = f.names
    if occursin("mbar", string(pcol))
        pressure_Pa = 100*f[pcol]
    else
        pressure_Pa = f[pcol]
    end
    if occursin("m3/h", string(scol))
        speed_m3ps = f[scol]/3600
    elseif occursin("l/s", string(scol))
        speed_m3ps = f[scol]/1000
    end
    spl = BSpline(log10.(pressure_Pa), speed_m3ps)
    p_Pa -> spl(log10.(p_Pa))
end

pumps = Dict(
    "nXDS15i" => get_characteristic("nXDS15i_m3ph_vs_mbar.csv"),
    "HiPace700" => get_characteristic("HiPace700_lps_vs_mbar.csv"),
    "HiPace80Neo" => get_characteristic("HiPace80Neo_lps_vs_mbar.csv")
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
        tube_PVflow(a, flength, gas, Pmax, p2) - p2*pump_speed(p2)
    end
end

Pascal_to_mbar(P_Pascal) = P_Pascal/100
slpm(qpv) = qpv/atm * 1000 * 60
mbarlps(qpv) = qpv / bar * 1000 * 1000

function choked_pressure(a, flength, gas, Pmax)
    ID = 2a
    η = gas_viscosity(gas)
    c = speed_of_sound(gas)
    ID^2*Pmax^2/(64η*flength*c)
end

function tube_PVflow_choked(a, flength, gas, Pmax)
    c = speed_of_sound(gas)
    pd = choked_pressure(a, flength, gas, Pmax)
    pd * π*a^2 * c
end

function gradient_end_pressure_choked(a, flength, gas, Pmax, pump_speed)
    qpv = tube_PVflow_choked(a, flength, gas, Pmax)
    find_zero(1e-3) do p2
        qpv - p2*pump_speed(p2)
    end
end

function gradient_end_pressure_choked(a, flength, gas, Pmax, pump_speed::Number)
    gradient_end_pressure_choked(a, flength, gas, Pmax, p -> pump_speed)
end

end