var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [HistogramsFit]\nOrder   = [:type, :function]","category":"page"},{"location":"api/#HistogramsFit.MultinomialBinsModel","page":"API","title":"HistogramsFit.MultinomialBinsModel","text":"MultinomialBinsModel{N,T,U,J} <: AbstractHistogramFitModel{N,J}\n\nThis struct stores all the informations needed to construct a χ² statistics where the bins respect a multinomial distribution.\n\nnote: Note\nThis statistics should be used whenever the total number of events\nis fixed by the measurement you performed.\n\n\n\n\n\n","category":"type"},{"location":"api/#HistogramsFit.PoissonianBinsModel","page":"API","title":"HistogramsFit.PoissonianBinsModel","text":"PoissonianBinsModel{N,T,U,J} <: AbstractHistogramFitModel{N,J}\n\nThis struct stores all the informations needed to construct a χ² statistics where all the bins have an independent poissonian distribution.\n\nnote: Note\nThis statistics should be used whenever the total number of events\nis **not** fixed by the measurement you performed.\n\n\n\n\n\n","category":"type"},{"location":"api/#HistogramsFit.hcubature-Tuple{Any, Any, Any}","page":"API","title":"HistogramsFit.hcubature","text":"    hcubature(f, domain, α)\n\nComputes n-d integrals using the HCubatureJL() method from Integrals.\n\n\n\n\n\n","category":"method"},{"location":"api/#HistogramsFit.quadgk-Tuple{Any, Any, Any}","page":"API","title":"HistogramsFit.quadgk","text":"    quadgk(f, domain, α)\n\nComputes 1-d integrals using the GuadGKJL() method from Integrals.\n\n\n\n\n\n","category":"method"},{"location":"api/#HistogramsFit.simpson-Tuple{Any, Any, Any}","page":"API","title":"HistogramsFit.simpson","text":"    simpson(f, domain, α)\n\nImplements the Simpson's rule to compute 1-d integrals numerically.\n\nIts usage is suggested only when the bins are narrow.\n\ntodo: Todo\nExplain how narrow is defined (compared to what?)\n\n\n\n\n\n","category":"method"},{"location":"#Welcome-to-HistrogramsFit.jl's-documentation","page":"Welcome to HistrogramsFit.jl's documentation","title":"Welcome to HistrogramsFit.jl's documentation","text":"","category":"section"},{"location":"","page":"Welcome to HistrogramsFit.jl's documentation","title":"Welcome to HistrogramsFit.jl's documentation","text":"The problem of fitting curves to histograms is ubiquitous in High Energy Physics (HEP) and usually it involves three steps:","category":"page"},{"location":"","page":"Welcome to HistrogramsFit.jl's documentation","title":"Welcome to HistrogramsFit.jl's documentation","text":"Determining the \"best fit\" parameters of a curve;\nDetermining the errors on the parameters;\nJudging the goodness of the fit. ","category":"page"},{"location":"","page":"Welcome to HistrogramsFit.jl's documentation","title":"Welcome to HistrogramsFit.jl's documentation","text":"This simple Julia module takes an nD histogram and a data distribution model and creates a theoretically sound chi square statistics that can be used to perform","category":"page"},{"location":"","page":"Welcome to HistrogramsFit.jl's documentation","title":"Welcome to HistrogramsFit.jl's documentation","text":"Point estimation;\nConfidence interval estimation;\nGoodness-of-fit testing.","category":"page"},{"location":"","page":"Welcome to HistrogramsFit.jl's documentation","title":"Welcome to HistrogramsFit.jl's documentation","text":"warning: Warning\nThis is still a work in progress and it is in pre-alpha stage. Expect issues and bug. If you want to contribute, please feel free to open a issue/start a discussion on GitHub! Thank you!","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"todo: Todo\nwrite the tutorials!!!","category":"page"}]
}
