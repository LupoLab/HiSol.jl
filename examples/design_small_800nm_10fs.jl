import HiSol

λ_target = 200e-9
gas = :HeJ

λ0 = 800e-9
τfwhm = 10e-15
energy = 100e-6

maxlength = 2

window = HiSol.WindowConstraint(λ0, :SiO2; thickness=1e-3, LIDT=2000)

fig, paramsf = HiSol.Design.aplot_energy_maxlength(λ_target, gas, λ0, τfwhm, energy, maxlength;
                                                   input_constraint=window,
                                                   output_constraint=window)
