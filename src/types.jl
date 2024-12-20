abstract type AbstractHistogramFitModel{N,J} end


"""
	PoissonianBinsModel{N,T,U,J} <: AbstractHistogramFitModel{N,J}

This struct stores all the informations needed to construct a
χ² statistics where all the bins have an independent poissonian
distribution.

!!! note
    This statistics should be used whenever the total number of events
    is **not** fixed by the measurement you performed. 

"""
@kwdef struct PoissonianBinsModel{N,T,U,J} <: AbstractHistogramFitModel{N,J}
        edges::NTuple{N,T}
        bincounts::AbstractArray{U,N}
        curve::SciMLBase.AbstractIntegralFunction
        params_names::NTuple{J,Symbol}
        integrator::SciMLBase.AbstractIntegralAlgorithm
end


"""
	MultinomialBinsModel{N,T,U,J} <: AbstractHistogramFitModel{N,J}

This struct stores all the informations needed to construct a
χ² statistics where the bins respect a multinomial distribution.

!!! note
    This statistics should be used whenever the total number of events
    is fixed by the measurement you performed. 

"""
@kwdef struct MultinomialBinsModel{N,T,U,J} <: AbstractHistogramFitModel{N,J}
        edges::NTuple{N,T}
        bincounts::AbstractArray{U,N}
        curve::SciMLBase.AbstractIntegralFunction
        params_names::NTuple{J,Symbol}
        integrator::SciMLBase.AbstractIntegralAlgorithm
end


for M in [:PoissonianBinsModel, :MultinomialBinsModel]
        @eval begin
                function $M(h::Histogram, fun, params_names; integrator=:default)
                        e = h.edges
                        n = h.weights
                        f = typeof(fun) <: SciMLBase.AbstractIntegralFunction ? fun : IntegralFunction(fun)
                        i = integrator
                        if i == :default
                                i = HCubatureJL()
                        end
                        return $M(edges=e, bincounts=n, curve=f, params_names=params_names, integrator=i)
                end
        end
        for H in [:Hist1D, :Hist2D, :Hist3D]
                @eval begin
                        function $M(h::$H, fun, params_names; integrator=:default)
                                e = $H != Hist1D ? Tuple(binedges(h)) : ([binedges(h)...],)
                                n = FHist.bincounts(h)
                                f = typeof(fun) <: SciMLBase.AbstractIntegralFunction ? fun : IntegralFunction(fun)
                                i = integrator
                                if i == :default
                                        i = HCubatureJL()
                                end
                                return $M(edges=e, bincounts=n, curve=f, params_names=params_names, integrator=i)
                        end
                end
        end
end
