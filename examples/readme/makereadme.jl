# To be run from the root of the repository to generate the README.md file.
using Literate
Literate.markdown("examples/readme/README.jl", "."; flavor = Literate.CommonMarkFlavor(), execute = true)