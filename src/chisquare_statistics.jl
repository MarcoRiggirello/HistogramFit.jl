abstract type AbstractChiSquareStatistics end


struct PoissonLikelihoodChiSquare{H,T,U} <: AbstractChiSquareStatistics
        data::H
        model::T
        χ²_λ::U
        α_keys
end


struct MultinomialLikelihoodChiSquare{H,T,U} <: AbstractChiSquareStatistics
        data::H
        model::T
        χ²_λ::U
        α_keys
end
