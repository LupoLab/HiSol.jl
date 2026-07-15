import CairoMakie #hide
using HiSol; #hide
wc = WindowConstraint(800e-9, :SiO2; thickness=1e-3, LIDT=2000); #hide
figs, _ = design_space_a_energy(200e-9, :He, 800e-9, 10e-15, 10; input_constraint=wc, output_constraint=wc); #hide

# > [!WARNING]
# > This package is a work in progress. Function signatures and internals are subject to change without notice. Please use with caution and contribute corrections or improvements if possible.

# > [!WARNING]
# > The way space constraints on the HCF system are specified has recently changed. See [API changes](#api-changes) for details and how to update existing scripts.

#=
# HiSol.jl

**HiSol.jl** is a Julia package to aid in the design of gas-filled hollow-core fibre systems for pulse compression and soliton dynamics using ultrashort laser pulses. It does this by implementing a collection of design rules developed in the [Laboratory of Ultrafast Physics and Optics](https://lupo-lab.com).

HiSol.jl is based on analytical expressions for the various parameters involved in soliton dynamics (e.g. the phase-matching wavelength for RDW emission). Since soliton self-compression is a highly dynamical process, any analytical model is only approximate. Design parameters found with HiSol.jl should be checked with full numerical simulations, e.g. using [Luna.jl](https://github.com/LupoLab/Luna.jl), before basing experimental designs on them.

## Installation
HiSol.jl is not yet a registered Julia package. You therefore need to install it directly from GitHub, using the following commands in the Julia REPL:
```julia
] add https://github.com/LupoLab/HiSol.jl
```
**or**
```julia
using Pkg
Pkg.add(url="https://github.com/LupoLab/HiSol.jl")
```
=#

#=
## Basic usage
The main functionality of HiSol.jl is to find the parameter space for soliton self-compression and especially resonant dispersive wave (RDW) emission in gas-filled hollow capillary fibres. The main function which shows the available design space is `design_space_a_energy`. It requires 5 input arguments and 2 keyword arguments:
```julia
design_space_a_energy(λ_target, gas, λ0, τfwhm, maxlength; input_constraint, output_constraint)
```
- `λ_target::Number`: the target wavelength for RDW emission in metres.
- `gas::Symbol`: the gas species.
- `λ0::Number`: the pump/driving wavelength in metres.
- `τfwhm::Number`: the pump/driving pulse duration in seconds
- `maxlength::Number`: the maximum **total** length of the HCF system in metres.
- `input_constraint`/`output_constraint`: length constraints for the entrance and exit side of the HCF (see below).

Because capillary fibres need to be kept perfectly straight, the maximum HCF length is determined by the available straight length of optical table. In most cases, space is required on both sides of the HCF to allow the incoming/outgoing beam to converge/diverge without being detrimentally affected by nonlinearities (in windows) or damaging the steering and focusing optics. How much space is set aside on each side of the HCF is determined by the length constraints (see [Length constraints](#length-constraints) below).

To create the plots, a plotting backend has to be loaded alongside HiSol.jl itself (see [Plotting backends](#plotting-backends) below). Here we use CairoMakie, which produces static plots; for interactive plot windows, load GLMakie instead.

As an example, we will design the HCF system used in the first demonstration of RDW emission in a hollow capillary fibre [Travers et al., Nature Photonics 13, 547 (2019)]. First we need to load the package and define our fixed parameters and constraints, then we call the function. Here we use a `WindowConstraint` for both ends of the HCF: a 1-mm thick fused-silica window which also has a damage threshold of 2000 J/m².
=#
using CairoMakie
using HiSol
dir = joinpath(pkgdir(HiSol), "examples/readme/figures/") #hide

λ_target = 160e-9 # 160 nm RDW
gas = :He # helium gas
λ0 = 800e-9 # 800 nm driving pulse
τfwhm = 10e-15 # 10 fs driving pulse
maxlength = 5 # 5 m maximum setup length

window = WindowConstraint(λ0, :SiO2; thickness=1e-3, LIDT=2000)

figs, params, as, energies, ratios = design_space_a_energy(λ_target, gas, λ0, τfwhm, maxlength;
                                                           input_constraint=window,
                                                           output_constraint=window)
save(joinpath(dir, "readme_ex_1a.svg"), figs[1]) #hide
save(joinpath(dir, "readme_ex_1b.svg"), figs[2]) #hide

# This will produce the following plots (the `Figure` objects are returned in the `figs` variable above.)
# ![Criteria ratios for design space example](examples/readme/figures/readme_ex_1a.svg)
# ![Fission length and soliton order in design space example](examples/readme/figures/readme_ex_1b.svg)

#=
The top row of plots shows parameter ratios which correspond to the four criteria we need to fulfill to observe RDW emission
- Loss: the fission length needs to be shorter than the $1/e$ loss length of the capillary. Equivalently, we need $\frac{L_\mathrm{fiss}}{L_\mathrm{loss}} < 1$.
- Fission length: the fission length needs to be shorter than the maximum capillary length which can fit into the available space for these parameters: $\frac{L_\mathrm{fiss}}{L_{\mathrm{HCF}}} < 1$.
- Minimum soliton order: for any soliton self-compression to occur, we need $N \geq 1.5$. Additionally, for RDWs far away from the pump wavelength (i.e. very short RDW wavelengths), we need $N \geq N_\mathrm{min}$ where the minimum soliton order $N_\mathrm{min}$ is based on an empirical relation. The relevant ratio is thus $\frac{N_\mathrm{min}}{N} < 1$.
- Maximum soliton order: we must not exceed the maximum soliton order, which encodes intensity limits due to photoionisation and plasma or self-focusing: $\frac{N}{N_\mathrm{max}} < 1$.

In each plot, the grey line shows the boundary between regions where the respective ratio is above and below $1$. The first plot also shows the resulting design space: the region where all four ratios are below $1$, and thus RDW emission should be possible.

The second figure shows the two key parameters for the soliton dynamics&mdash;the fission length and the soliton order&mdash;within the design space outlined by the four criteria. 

To be more precise in our choices, instead of reading numbers off of the plot, we can use the `params` function returned by `design_space_a_energy`. This takes the two coordinates of the figure (core radius, pulse energy) and returns a full list of the specifications of the system. For example, we can find the exact configuration for a 125 μm core radius and 200 μJ:
=#
a = 125e-6 # 125 μm
energy = 200e-6 # 200 μJ
p = params(a, energy)
#=
Here `p` is now a `NamedTuple` containing the parameters. Its fields are:
- `radius`: (input) the chosen core radius.
- `energy`: (input) the chosen pulse energy.
- `density`: gas density (m⁻³).
- `pressure`: pressure (bar).
- `Lfiss`: fission length.
- `N`: soliton order.
- `flength`: *maximum* fibre length which fits into the space for the chosen parameters.
- `τfwhm`: FWHM pulse duration.
- `Nmin`: minimum soliton order.
- `Nmax`: maximum soliton order.
- `Lloss`: loss length of the HCF.
- `Isupp`: barrier suppression intensity of the filling gas.
- `Icrit`: critical intensity for self-focusing, calculated as $I_\mathrm{crit} = P_\mathrm{crit}/A_\mathrm{eff}$, where $P_\mathrm{crit}$ is the critical *power* and $A_\mathrm{eff}$ is the effective area of the waveguide.
- `intensity`: peak intensity of the driving pulse (W/m²), calculated as $P_0/A_\mathrm{eff}$ where $P_0$ is the pulse peak power.

To find, for example, the gas pressure for one point in the design space, we can access the respective field:
=#
p.pressure
# Or, similarly, the fission length:
p.Lfiss

#=
## Additional options
The `design_space_a_energy` function takes a large range of additional keyword arguments which affect the design rules it applies. The two main categories concern a) safety factors for the nonlinear dynamics themselves b) limitations on the maximum HCF length in the given space (`maxlength`).

### Safety factors
The three main keyword arguments which affect the "safety margin" in the calculations control safety factors. For each safety factor, a *larger* number implies a *more conservative* rule, which implies a smaller design space.
- `S_sf`: safety factor on the critical power for self-focusing. The driving-pulse peak power will not exceed $P_\mathrm{crit}/S_\mathrm{sf}$. (Default: 5)
- `S_ion`: safety factor on strong-field photoionisation. The driving-pulse intensity (in the mode-averaged sense, i.e. $P_0/A_\mathrm{eff}$) will not exceed $I_\mathrm{supp}/S_\mathrm{ion}$, where $I_\mathrm{supp}$ is the barrier suppression intensity of the gas. (Default: 10)
- `S_fiss`: safety factor on the fission length. The *minimum* fibre length for efficient RDW emission is taken to be $S_\mathrm{fiss}L_\mathrm{fiss}$, where $L_\mathrm{fiss}$ is the fission length. (Default: 1.5)

`S_fiss` affects both the loss cut-off and the length cut-off: for a parameter combination to be included in the design space, $S_\mathrm{fiss}\times L_\mathrm{fiss}$ needs to be shorter than both the maximum HCF length *and* the loss length. A fourth keyword argument can alter that behaviour:

- `S_loss`: safety factor on the loss length. The fission length (multiplied by the safety factor) will be shorter than $L_\mathrm{loss}/S_\mathrm{loss}$. (Default: 1)

Note that `S_loss` follows the same convention as the other safety factors: a larger value results in a more conservative design rule.

The default safety factors are intentionally very conservative, with the aim of generating parameter combinations which are very likely to work in practice. In extreme cases, it can be necessary to adjust the safety factors. For example, RDW emission at very short wavelengths is commonly ionisation-limited. To achieve efficient frequency conversion, the conservative limit can be exceeded by several times&mdash;at the cost of relying on more extreme nonlinear dynamics which can be more sensitive to minor perturbations. Similarly, the loss limit can be exceeded (see e.g. Chen *et al.*, 10.1364/OL.553345) at the cost of a potentially significant reduction in conversion efficiency to the RDW. 

### Length constraints
How much of the available space is taken up by the free-space propagation at either end of the HCF is defined by the required keyword arguments `input_constraint` and `output_constraint`. Each of these is a callable object—a subtype of `HiSol.Focusing.LengthConstraint`—which returns the minimum distance between the end of the HCF and the first/last other optical element. Because the two ends are specified separately, they can have completely different constraints (e.g. different window materials, or a window at the entrance but an open exit towards a mirror). Four types of constraint are available:

- `NoConstraint()`: no space is required at all—the HCF can take up the whole available length.
- `FixedConstraint(distance)`: a fixed `distance` is set aside, regardless of all other parameters.
- `DamageConstraint(λref, LIDT; S_fluence=5, conversion=1)`: keeps the fluence on the first/last mirror below the damage threshold `LIDT` (**in SI units**, i.e. J/m²) divided by the safety factor `S_fluence`. `conversion` is an energy conversion factor between the pulse energy in the HCF and the energy hitting the mirror (e.g. a frequency-conversion efficiency on the exit side). Some example values for `LIDT`:
    - Newport ultrafast silver: 0.12 J/cm² = 1200 J/m²
    - Thorlabs protected silver: 0.225 J/cm² = 2250 J/m²
    - Layertec fs-optimised protected silver: 0.38 J/cm² = 3800 J/m²
- `WindowConstraint(λref, n2; kwargs...)`: a window made of a material with nonlinear refractive index `n2` (either a `Symbol` like `:SiO2` or `:MgF2`, or a `Number` in m²/W). The constraint keeps the B-integral in the window below `Bmax` (default 0.2, see 10.1364/OE.482749), while requiring that the window is thick enough to hold the gas pressure and large enough to pass the beam. The window `thickness` and `aperture` radius can each be given as fixed numbers or left variable, in which case suitable values are found automatically. A window damage threshold can also be taken into account by passing `LIDT`. See the docstring of `WindowConstraint` for the full list of keyword arguments.

The distance, thickness and aperture implied by a constraint for a given HCF design can be inspected with `HiSol.Focusing.details(constraint, a, energy, τfwhm; pressure)`.

Again, conservative defaults (e.g. the safety factors `S_fluence` and `S_break`) are used with the aim of producing working parameter combinations.

## Plotting backends
HiSol.jl does not itself depend on a plotting package. To create plots, load one of the supported backends alongside HiSol.jl:

- [GLMakie](https://docs.makie.org/): interactive plot windows and the **interactive design-space explorers** (see below). Recommended for interactive use.
- [CairoMakie](https://docs.makie.org/): static publication-quality output (SVG/PDF/PNG).
- [PythonPlot](https://github.com/JuliaPy/PythonPlot.jl) or [PyPlot](https://github.com/JuliaPy/PyPlot.jl): the full capability of matplotlib.

For example:
```julia
using GLMakie
using HiSol
# ... plotting functions are now available
```
Calling a plotting function without a backend loaded throws an error. The interactive explorers and the 3D design-space plot (below) are only available with Makie backends.

## Interactive and 3D design-space exploration
With GLMakie loaded, `interactive_design_space()` opens a window showing the same maps as `design_space_a_energy`, with graphical controls for all parameters (target wavelength, gas, pump wavelength and duration, setup length, safety factors and the length constraints at both ends). Press *Compute* to regenerate the maps and click in any map to inspect the full system specification at that point.

The design space can also be shown in three dimensions, with core radius and energy as two of the axes and either the target wavelength or the pulse duration as the third. Passing a range for exactly one of `λ_target` or `τfwhm` selects the third axis:
```julia
fig, data = design_space_a_energy_3D(range(120e-9, 300e-9, 32), gas, λ0, τfwhm, maxlength;
                                     input_constraint=window, output_constraint=window)
```
The boundary of the design space is drawn as an isosurface, which can be rotated freely in GLMakie. The interactive version, `interactive_design_space_3D()`, additionally shows a 2D slice which can be moved along the third axis with a slider. See [`examples/design_space_gui.jl`](examples/design_space_gui.jl), [`examples/design_space_3D.jl`](examples/design_space_3D.jl) and [`examples/design_space_3D_gui.jl`](examples/design_space_3D_gui.jl).

## API changes
The way space constraints on the HCF system are specified has recently changed. Previously, `design_space_a_energy` and the functions in `HiSol.Compressor` took a fixed set of keyword arguments (`entrance_window`, `exit_window`, `thickness`, `material`, `Bmax`, `LIDT` and `S_fluence`) which described one window geometry and one mirror damage threshold for both ends of the HCF. These keyword arguments have been **removed** and replaced by the two required keyword arguments `input_constraint` and `output_constraint`, each of which can be any of the constraint types described under [Length constraints](#length-constraints). This makes it possible to specify completely different constraints for the entrance and exit sides, to let window thickness and aperture be found automatically, and to add custom constraint types.

To update existing scripts:

- The old **defaults** (`entrance_window=true, exit_window=true, thickness=1e-3, material=:SiO2, Bmax=0.2, LIDT=2000, S_fluence=5`) are equivalent to passing
  ```julia
  window = WindowConstraint(λ0, :SiO2; thickness=1e-3, Bmax=0.2, LIDT=2000, S_fluence=5)
  design_space_a_energy(λ_target, gas, λ0, τfwhm, maxlength;
                        input_constraint=window, output_constraint=window)
  ```
  Passing `LIDT` to the `WindowConstraint` reproduces the old behaviour, in which the mirror damage threshold was taken into account even when a window was present.
- `entrance_window=false`/`exit_window=false` (no window, mirror damage threshold only) is equivalent to using `DamageConstraint(λ0, LIDT; S_fluence)` for the respective side.

Other changes:

- The function `params` returned by `design_space_a_energy` still takes `(a, energy)` as before, and now optionally accepts the pulse duration as a third argument.
- `HiSol.Focusing.max_flength` is now called as `max_flength(a, energy, τfwhm, maxlength, input_constraint, output_constraint; pressure)`.
- `HiSol.Compressor.optimise` and `HiSol.Compressor.plot_optimise` take the same new keyword arguments. The parameters returned by `HiSol.Compressor.params_maxlength` contain the entrance/exit distances `d_in`/`d_out` and beam sizes `beamsize_w0_in`/`beamsize_w0_out` instead of the symmetric `window_distance` and `beamsize_w0`.
- `HiSol.Design.max_energy`, `HiSol.Design.maximum_radius`, and `HiSol.Design.aplot_energy_maxlength` also take `input_constraint`/`output_constraint` instead of the old keyword arguments.
- `HiSol.Design.params_maxlength` has been removed; its role is taken over by `HiSol.Design.Params`.
- HiSol.jl no longer depends on PyPlot directly. A plotting backend (GLMakie, CairoMakie, PythonPlot or PyPlot) must be loaded to create plots—see [Plotting backends](#plotting-backends).
=#