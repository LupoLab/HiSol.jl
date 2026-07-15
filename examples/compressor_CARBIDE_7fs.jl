using GLMakie # load a plotting backend (or use PythonPlot / CairoMakie)
import HiSol

λ0 = 1030e-9
τfwhm_in = 29e-15
τfwhm_out = 7e-15

energyout = 511e-6
a0 = 225e-6
flength0 = 1.8
tr0 = HiSol.HCF.transmission(flength0, a0, λ0)

energy = 650e-6

gas = :Ar
maxlength = 3.56

window = HiSol.WindowConstraint(λ0, :MgF2; thickness=1e-3)

HiSol.Compressor.plot_optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                               input_constraint=window, output_constraint=window,
                               S_sf=2.5, S_ion=10)
