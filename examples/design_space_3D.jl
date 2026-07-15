# Three-dimensional representation of the design space (requires a Makie backend;
# use GLMakie for a rotatable interactive view). The design-space boundary is drawn as
# an isosurface over core radius, pulse energy, and a third axis, which is chosen by
# passing a range for either the target wavelength or the pulse duration.
using GLMakie
using HiSol

gas = :He
λ0 = 800e-9
τfwhm = 10e-15
maxlength = 5

window = WindowConstraint(λ0, :SiO2; thickness=1e-3, LIDT=2000)

# third axis: RDW target wavelength from 120 nm to 300 nm
fig1, data1 = design_space_a_energy_3D(range(120e-9, 300e-9, 32), gas, λ0, τfwhm, maxlength;
                                       input_constraint=window, output_constraint=window)

# third axis: pump pulse duration from 5 fs to 25 fs at a fixed 160 nm target
fig2, data2 = design_space_a_energy_3D(160e-9, gas, λ0, range(5e-15, 25e-15, 32), maxlength;
                                       input_constraint=window, output_constraint=window)
