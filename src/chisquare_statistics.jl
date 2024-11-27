function predicted_events(hfm::AbstractHistogramFitModel{N,J}, bin::CartesianIndex{N}) where {N,J}
        f = curve(hfm)
        q = integrator(hfm)
        e = edges(hfm)
        l = getindex.(e, bin)
        u = getindex.(e, bin + one(bin))
        return q(f, l, u)
end


function bin_lnλ(hfm::PoissonianBinsModel{N,T,U,J}, bin::CartesianIndex{N}) where {N,T,U,J}
        n = bincounts(hfm)[bin]
        y = predicted_events(hfm, bin)
        return y - n + n * log(n / y)
end


function bin_lnλ(hfm::MultinomialBinsModel{N,T,U,J}, bin::CartesianIndex{N}) where {N,T,U,J}
        n = bincounts(hfm)[bin]
        y = predicted_events(hfm, bin)
        return n * log(n / y)
end


function chisquare(hfm::AbstractHistogramFitModel)
        χ² = 0
        for b in CartesianIndices(bincounts(hfm))
                χ² += bin_lnλ(hfm, b)
        end
        return 2χ²
end
