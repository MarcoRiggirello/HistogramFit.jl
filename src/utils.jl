curve(hfm::AbstractHistogramFitModel) = curve(hfm)

curve(hfm::PoissonianBinsModel) = hfm.curve
curve(hfm::MultinomialBinsModel) = hfm.curve


integrator(hfm::AbstractHistogramFitModel) = integrator(hfm)

integrator(hfm::PoissonianBinsModel) = hfm.integrator
integrator(hfm::MultinomialBinsModel) = hfm.integrator


bincounts(hfm::AbstractHistogramFitModel) = bincounts(hfm)

bincounts(hfm::PoissonianBinsModel) = hfm.bincounts
bincounts(hfm::MultinomialBinsModel) = hfm.bincounts


edges(hfm::AbstractHistogramFitModel) = edges(hfm)

edges(hfm::PoissonianBinsModel) = hfm.edges
edges(hfm::MultinomialBinsModel) = hfm.edges
