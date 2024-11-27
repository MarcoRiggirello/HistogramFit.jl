module HistogramsFit

include("types.jl")
include("utils.jl")
include("simpson.jl")
include("chisquare_statistics.jl")

export PoissonianBinsModel, MultinomialBinsModel
export chisquare
export simpson

end
