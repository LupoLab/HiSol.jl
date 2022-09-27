using HiSol
using Documenter

DocMeta.setdocmeta!(HiSol, :DocTestSetup, :(using HiSol); recursive=true)

makedocs(;
    modules=[HiSol],
    authors="Christian Brahms <c.brahms@hw.ac.uk>",
    repo="https://github.com/LupoLab/HiSol.jl/blob/{commit}{path}#{line}",
    sitename="HiSol.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://LupoLab.github.io/HiSol.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LupoLab/HiSol.jl",
    devbranch="master",
)
