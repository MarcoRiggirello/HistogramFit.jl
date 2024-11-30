push!(LOAD_PATH,"../src/")

using Documenter, HistogramsFit

makedocs(
        sitename = "HistogramsFit.jl",
        modules = [HistogramsFit],
	pages = ["index.md",
		 "tutorial.md",
		 "api.md"]
)

deploydocs(
    repo = "github.com/MarcoRiggirello/HistogramsFit.jl.git",
)
