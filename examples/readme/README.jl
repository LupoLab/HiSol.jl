import PyPlot: plt, pygui #hide
pygui(false); #hide

# > [!WARNING]
# > This package is a work in progress. Function signatures and internals are subject to change without notice. Please use with caution and contribute corrections or improvements if possible. 

#=
# HiSol.jl

**HiSol.jl** is a Julia package to aid in the design of gas-filled hollow-core fibre systems for pulse compression and soliton dynamics using ultrashort laser pulses. It does this by implementing a collection of design rules developed in the [Laboratory of Ultrafast Physics and Optics](https://lupo-lab.com).

## Installation
HiSol.jl is not yet a registered Julia package. You therefore need to nstall it directly from GitHub, using the following commands in the Julia REPL:
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
The main functionality of HiSol.jl is to find the parameter space for soliton self-compression and especially resonant dispersive wave (RDW) emission in gas-filled hollow capillary fibres. The main function which shows the available design space is `design_space_a_energy`. It requires 5 input arguments:
```julia
design_space_a_energy(λ_target, gas, λ0, τfwhm, maxlength)
```
- `λ_target::Number`: the target wavelength for RDW emission in metres.
- `gas::Symbol`: the gas species.
- `λ0::Number`: the pump/driving wavelength in metres.
- `τfwhm::Number`: the pump/driving pulse duration in seconds
- `maxlength::Number`: the maximum **total** length of the HCF system in metres.

As an example, we will design the HCF system used in the first demonstration of RDW emission in a hollow capillary fibre [Travers et al., Nature Photonics 13, 547 (2019)]. First we need to load the package and define our fixed parameters and constraints, then we call the function.
=#
using HiSol;
dir = joinpath(pkgdir(HiSol), "examples/readme/figures/") #hide

λ_target = 160e-9 # 160 nm RDW
gas = :He # helium gas
λ0 = 800e-9 # 800 nm driving pulse 
τfwhm = 10e-15 # 10 fs driving pulse
maxlength = 5 # 5 m maximum setup length

figs, params, a, energy, ratios = design_space_a_energy(λ_target, gas, λ0, τfwhm, maxlength)
figs[1].savefig(joinpath(dir, "readme_ex_1a.svg"))
figs[2].savefig(joinpath(dir, "readme_ex_1b.svg"))

# This will produce the following plots (the `Figure` objects returned in the `figs` variable above.)
# ![Criteria ratios for design space example](examples/readme/figures/readme_ex_1a.svg)
# ![Fission length and soliton order in design space example](examples/readme/figures/readme_ex_1b.svg)