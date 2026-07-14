module HiSol

include("Data.jl")
include("HCF.jl")
include("Solitons.jl")
include("Limits.jl")
include("Focusing.jl")
include("Compressor.jl")
include("GasFlow.jl")
include("Design.jl")

import .Design: design_space_a_energy
export design_space_a_energy

end
