import HiSol

λ0 = 1030e-9
τfwhm_in = 220e-15
τfwhm_out = 50e-15

energy = 200e-6 # 100 kHz, 80 W

gas = :Kr
maxlength = 1.5

window = HiSol.WindowConstraint(λ0, :MgF2; thickness=1e-3)

HiSol.Compressor.plot_optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                               input_constraint=window, output_constraint=window)
