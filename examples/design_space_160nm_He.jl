# Design-space example from the README: the HCF system used in the first demonstration
# of RDW emission in a hollow capillary fibre [Travers et al., Nature Photonics 13, 547 (2019)].
using GLMakie # load a plotting backend (or use PythonPlot / CairoMakie)
using HiSol

λ_target = 160e-9 # 160 nm RDW
gas = :He # helium gas
λ0 = 800e-9 # 800 nm driving pulse
τfwhm = 10e-15 # 10 fs driving pulse
maxlength = 5 # 5 m maximum setup length

# 1-mm thick fused-silica window with a damage threshold of 2000 J/m² at both ends
window = WindowConstraint(λ0, :SiO2; thickness=1e-3, LIDT=2000)

figs, params, as, energies, ratios = design_space_a_energy(λ_target, gas, λ0, τfwhm, maxlength;
                                                           input_constraint=window,
                                                           output_constraint=window)

# full system specification at one point in the design space
a = 125e-6 # 125 μm core radius
energy = 200e-6 # 200 μJ
p = params(a, energy)
display(p)
