import HISOL

λ_target = 200e-9
gas = :HeJ

λ0 = 800e-9
τfwhm = 10e-15
energy = 50e-6

maxlength = 2

fig, paramsf = HISOL.Design.aplot_energy_maxlength(λ_target, gas, λ0, τfwhm, energy, maxlength)
