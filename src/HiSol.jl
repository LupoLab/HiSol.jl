module HiSol

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

import .Design: design_space_a_energy
export design_space_a_energy

end
