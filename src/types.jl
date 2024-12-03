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
        curve::Function
        params_names::NTuple{J,Symbol}
        integrator::Function
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
        curve::Function
        params_names::NTuple{J,Symbol}
        integrator::Function
end


for M in [:PoissonianBinsModel, :MultinomialBinsModel]
        @eval begin
                function $M(h::Histogram, f, params_names; integrator=:default)
                        e = h.edges
                        n = h.weights
                        i = integrator
                        if i == :default
                                i = length(e) > 1 ? hcubature : quadgk
                        end
                        return $M(edges=e, bincounts=n, curve=f, params_names=params_names, integrator=i)
                end
        end
        @eval begin
                function $M(h::Hist1D, f, params_names; integrator=:default)
                        e = ([binedges(h)...],)
                        n = FHist.bincounts(h)
                        i = integrator
                        if i == :default
                                i = quadgk
                        end
                        return $M(edges=e, bincounts=n, curve=f, params_names=params_names, integrator=i)
                end
        end
        for H in [:Hist2D, :Hist3D]
                @eval begin
                        function $M(h::$H, f, params_names; integrator=:default)
                                e = Tuple(binedges(h))
                                n = FHist.bincounts(h)
                                i = integrator
                                if i == :default
                                        i = hcubature
                                end
                                return $M(edges=e, bincounts=n, curve=f, params_names=params_names, integrator=i)
                        end
                end
        end
end
