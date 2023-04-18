using BioRecordsProcessing
using Documenter

DocMeta.setdocmeta!(BioRecordsProcessing, :DocTestSetup, :(using BioRecordsProcessing); recursive=true)

makedocs(;
    modules=[BioRecordsProcessing],
    authors="Jonathan Bieler <jonathan.bieler@alumni.epfl.ch> and contributors",
    repo="https://github.com/jonathanBieler/BioRecordsProcessing.jl/blob/{commit}{path}#{line}",
    sitename="BioRecordsProcessing.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jonathanBieler.github.io/BioRecordsProcessing.jl",
        assets=String[],
    ),
    pages=[
        "Manual" => "index.md",
        "Examples" => "examples.md",
        #"API" => "API.md",
    ],
)

deploydocs(;
    repo="github.com/jonathanBieler/BioRecordsProcessing.jl",
    devbranch = "main",
)
