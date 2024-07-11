using NeidSolarScripts
using Documenter

DocMeta.setdocmeta!(NeidSolarScripts, :DocTestSetup, :(using NeidSolarScripts); recursive=true)

makedocs(;
    modules=[NeidSolarScripts],
    authors="Eric Ford",
    repo="https://github.com/RvSpectML/NeidSolarScripts.jl/blob/{commit}{path}#{line}",
    sitename="NeidSolarScripts.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://RvSpectML.github.io/NeidSolarScripts.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
   #checkdocs=:none,
   checkdocs=:exports,
)

deploydocs(;
    repo="github.com/RvSpectML/NeidSolarScripts.jl",
    devbranch="main",
)
