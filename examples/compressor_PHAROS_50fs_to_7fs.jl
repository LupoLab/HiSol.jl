import HISOL

λ0 = 1030e-9
τfwhm_in = 50e-15
τfwhm_out = 7e-15

energy = 150e-6 # 100 kHz, 80 W

gas = :Ar
maxlength = 2.4

thickness = 1e-3
material = :MgF2

HISOL.Compressor.plot_optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength; S_sf=7)
 