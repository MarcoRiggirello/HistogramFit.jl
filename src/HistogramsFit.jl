module HistogramsFit

using StatsBase, FHist
using SciMLBase, Integrals

#include("integrators.jl")
include("types.jl")
include("utils.jl")
include("chisquare_statistics.jl")

export PoissonianBinsModel, MultinomialBinsModel
export chisquare_statistics

end
