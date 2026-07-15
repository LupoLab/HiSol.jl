# Interactive 3D design-space explorer. Requires an interactive Makie backend (GLMakie).
# Shows the design-space boundary as a rotatable isosurface plus a 2D slice which can be
# moved along the third axis with the slider. The third axis is chosen by passing a range
# for either the target wavelength (as here) or the pulse duration; the other parameters
# can be changed graphically. Press "Compute" to regenerate the volume.
using GLMakie
using HiSol

λ0 = 800e-9
window = WindowConstraint(λ0, :SiO2; thickness=1e-3, LIDT=2000)

interactive_design_space_3D(;
    λ_target=range(120e-9, 300e-9, 32), gas=:He, λ0, τfwhm=10e-15, maxlength=5,
    input_constraint=window, output_constraint=window)
