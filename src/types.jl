abstract type AbstractHistogramFitModel{N,J} end


@kwdef struct PoissonianBinsModel{N,T,U,J} <: AbstractHistogramFitModel{N,J}
        edges::NTuple{N,T}
        bincounts::AbstractArray{U,N}
        curve::Function
        params_names::NTuple{J,Symbol}
        integrator::Function
end


@kwdef struct MultinomialBinsModel{N,T,U,J} <: AbstractHistogramFitModel{N,J}
        edges::NTuple{N,T}
        bincounts::AbstractArray{U,N}
        curve::Function
        params_names::NTuple{J,Symbol}
        integrator::Function
end
