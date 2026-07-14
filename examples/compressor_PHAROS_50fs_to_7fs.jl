import HiSol

λ0 = 1030e-9
τfwhm_in = 50e-15
τfwhm_out = 7e-15

energy = 150e-6 # 100 kHz, 80 W

gas = :Ar
maxlength = 2.4

window = HiSol.WindowConstraint(λ0, :MgF2; thickness=1e-3)

HiSol.Compressor.plot_optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                               input_constraint=window, output_constraint=window, S_sf=7)
