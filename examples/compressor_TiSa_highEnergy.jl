import HiSol

λ0 = 800e-9
τfwhm_in = 35e-15
τfwhm_out = 10e-15

energy = 6e-3

gas = :He
maxlength = 8

# no windows; the limit is the damage threshold of the mirrors
mirror = HiSol.DamageConstraint(λ0, 3500; S_fluence=5)

HiSol.Compressor.plot_optimise(τfwhm_in, τfwhm_out, gas, λ0, energy, maxlength;
                               input_constraint=mirror, output_constraint=mirror, S_sf=5)
