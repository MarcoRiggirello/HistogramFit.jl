push!(LOAD_PATH, "../src/")

using Documenter, Literate
using HistogramsFit

# generate examples
TUTORIAL_DIR = joinpath(@__DIR__, "tutorials")
OUTDIR = joinpath(@__DIR__, "src", "tutorials_gen")
isdir(OUTDIR) && rm(OUTDIR, recursive=true)
mkpath(OUTDIR)

for tutorial in filter(contains(r".jl$"), readdir(TUTORIAL_DIR, join=true))
        Literate.markdown(tutorial, OUTDIR)
end

makedocs(
        sitename="HistogramsFit.jl",
        authors="Marco Riggirello",
        modules=[HistogramsFit],
        pages=["Home" => "index.md",
                "Tutorials" => Any[
                        "tutorials_gen/simple.md",
                        "tutorials_gen/hep.md",
                ],
                "Reference" => "api.md"
        ]
)

deploydocs(
        repo="github.com/MarcoRiggirello/HistogramsFit.jl.git",
)
