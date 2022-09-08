import HISOL

λ0 = 800e-9
τfwhm_in = 35e-15
τfwhm_out = 10e-15

energy = 6e-3

gas = :He
maxlength = 8

HISOL.Compressor.plot_optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                               entrance_window=false, exit_window=false, LIDT=3500, S_fluence=5, S_sf=5)
