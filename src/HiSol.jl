module HiSol

include("Plotting.jl")
include("Data.jl")
include("HCF.jl")
include("Solitons.jl")
include("Limits.jl")
include("Focusing.jl")
include("Compressor.jl")
include("GasFlow.jl")
include("Design.jl")

import .Focusing: LengthConstraint, NoConstraint, FixedConstraint, DamageConstraint, WindowConstraint
export NoConstraint, FixedConstraint, DamageConstraint, WindowConstraint

import .Design: design_space_a_energy, design_space_a_energy_3D,
                interactive_design_space, interactive_design_space_3D
export design_space_a_energy, design_space_a_energy_3D,
       interactive_design_space, interactive_design_space_3D

end
