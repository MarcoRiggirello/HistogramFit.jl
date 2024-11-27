abstract type AbstractHistogramFitModel{N<:Integer,J<:Integer} end


struct PoissonianBinsModel{N,T,U,J} <: AbstractHistogramFitModel{N,J}
        edges::NTuple{N,T}
        bincounts::AbstractArray{U,N}
        curve::Function
        params_names::NTuple{J,Symbol}
        integrator::Function
end


struct MultinomialBinsModel{N,T,U,J} <: AbstractHistogramFitModel{N,J}
        edges::NTuple{N,T}
        bincounts::AbstractArray{U,N}
        curve::Function
        params_names::NTuple{J,Symbol}
        integrator::Function
end
