using HISOL
using Documenter

DocMeta.setdocmeta!(HISOL, :DocTestSetup, :(using HISOL); recursive=true)

makedocs(;
    modules=[HISOL],
    authors="Christian Brahms <c.brahms@hw.ac.uk>",
    repo="https://github.com/LupoLab/HISOL.jl/blob/{commit}{path}#{line}",
    sitename="HISOL.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://LupoLab.github.io/HISOL.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LupoLab/HISOL.jl",
    devbranch="master",
)
