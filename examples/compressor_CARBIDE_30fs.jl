import HISOL

λ0 = 1030e-9
τfwhm_in = 330e-15
τfwhm_out = 30e-15

energy = 800e-6 # 100 kHz, 80 W

gas = :Ar
maxlength = 2.5

thickness = 2e-3
material = :MgF2

HISOL.Compressor.plot_optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength)
