abstract type AbstractHistogramFitModel{N,J} end

"""
	PoissonianBinsModel{N,T,U,J} <: AbstractHistogramFitModel{N,J}

This struct stores all the informations needed to construct a
χ² statistics where all the bins have an independent poissonian
distribution.

```
!!! note
        This statistics should be used whenever the total number of events
	is **not** fixed by the measurement you performed. 
```

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

```
!!! note
        This statistics should be used whenever the total number of events
	is fixed by the measurement you performed. 
```

"""
@kwdef struct MultinomialBinsModel{N,T,U,J} <: AbstractHistogramFitModel{N,J}
        edges::NTuple{N,T}
        bincounts::AbstractArray{U,N}
        curve::Function
        params_names::NTuple{J,Symbol}
        integrator::Function
end
