module HistogramsFit

using StatsBase, FHist
using Integrals

include("integrators.jl")
include("types.jl")
include("utils.jl")
include("chisquare_statistics.jl")

export PoissonianBinsModel, MultinomialBinsModel
export chisquare
export simpson

end
