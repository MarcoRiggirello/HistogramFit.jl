function predicted_events(α::AbstractVector, hfm::AbstractHistogramFitModel{N,J}, bin::CartesianIndex{N}) where {N,J}
        e = edges(hfm)
        l = [e[i][j] for (i, j) in enumerate(Tuple(bin))]
        u = [e[i][j] for (i, j) in enumerate(Tuple(bin + oneunit(bin)))]
        f = curve(hfm)
        q = integrator(hfm)
        p = IntegralProblem(f, (l, u), α)
        return solve(p, q).u
end


function bin_lnλ(α, hfm::PoissonianBinsModel{N,T,U,J}, bin::CartesianIndex{N}) where {N,T,U,J}
        n = bincounts(hfm)[bin]
        y = predicted_events(α, hfm, bin)
        return n == 0 ? y : y - n + n * log(n / y)
end


function bin_lnλ(α, hfm::MultinomialBinsModel{N,T,U,J}, bin::CartesianIndex{N}) where {N,T,U,J}
        n = bincounts(hfm)[bin]
        y = predicted_events(α, hfm, bin)
        return n == 0 ? 0 : n * log(n / y)
end


function chisquare_statistics(α, hfm::AbstractHistogramFitModel)
        lnλ = 0
        for b in CartesianIndices(bincounts(hfm))
                lnλ += bin_lnλ(α, hfm, b)
        end
        return 2lnλ
end
