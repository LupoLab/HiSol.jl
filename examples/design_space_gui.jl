# Interactive design-space explorer. Requires an interactive Makie backend (GLMakie).
# All parameters (target wavelength, gas, pump wavelength and duration, setup length,
# safety factors and the length constraints at both ends) can be changed graphically;
# press "Compute" to regenerate the maps. Click in any map to inspect a design point.
using GLMakie
using HiSol

λ0 = 800e-9
window = WindowConstraint(λ0, :SiO2; thickness=1e-3, LIDT=2000)

interactive_design_space(;
    λ_target=160e-9, gas=:He, λ0, τfwhm=10e-15, maxlength=5,
    input_constraint=window, output_constraint=window)
