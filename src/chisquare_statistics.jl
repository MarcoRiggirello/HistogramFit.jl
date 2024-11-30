function predicted_events(hfm::AbstractHistogramFitModel{N,J}, α::AbstractVector, bin::CartesianIndex{N}) where {N,J}
        f = curve(hfm)
        q = integrator(hfm)
        e = edges(hfm)
        l = [e[i][j] for (i, j) in enumerate(Tuple(bin))]
        u = [e[i][j] for (i, j) in enumerate(Tuple(bin + oneunit(bin)))]
        return q(f, (l, u), α)
end


function bin_lnλ(hfm::PoissonianBinsModel{N,T,U,J}, α, bin::CartesianIndex{N}) where {N,T,U,J}
        n = bincounts(hfm)[bin]
        y = predicted_events(hfm, α, bin)
        return y - n + n * log(n / y)
end


function bin_lnλ(hfm::MultinomialBinsModel{N,T,U,J}, α, bin::CartesianIndex{N}) where {N,T,U,J}
        n = bincounts(hfm)[bin]
        y = predicted_events(hfm, α, bin)
        return n * log(n / y)
end


function chisquare(hfm::AbstractHistogramFitModel, α)
        lnλ = 0
        for b in CartesianIndices(bincounts(hfm))
                lnλ += bin_lnλ(hfm, α, b)
        end
        return 2lnλ
end


χ²(hfm, α) = chisquare(hfm, α)
