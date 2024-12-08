function predicted_events(hfm::AbstractHistogramFitModel{N,J}, α::AbstractVector, bin::CartesianIndex{N}) where {N,J}
        e = edges(hfm)
        l = [e[i][j] for (i, j) in enumerate(Tuple(bin))]
        u = [e[i][j] for (i, j) in enumerate(Tuple(bin + oneunit(bin)))]
        f = curve(hfm)
        q = integrator(hfm)
        p = IntegralProblem(f, (l, u), α)
        return solve(p, q).u
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


function chisquare_statistics(hfm::AbstractHistogramFitModel, α)
        @floop for b in CartesianIndices(bincounts(hfm))
                blnλ = bin_lnλ(hfm, α, b)
                @reduce lnλ += blnλ
        end
        return 2lnλ
end
#function sum_bin_lnλ(hfm::AbstractHistogramFitModel, α, bins)
#        lnλ = 0
#        for b in bins
#                lnλ += bin_lnλ(hfm, α, b)
#        end
#        return lnλ
#end
#
#
#function chisquare_statistics(hfm::AbstractHistogramFitModel, α)
#        bins = CartesianIndices(bincounts(hfm))
#        T = typeof(bin_lnλ(hfm, α, bins[1]))
#        chunks = Iterators.partition(bins, length(bins) ÷ Threads.nthreads())
#        tasks = map(chunks) do chunk
#                Threads.@spawn sum_bin_lnλ(hfm, α, chunk)
#        end
#        lnλ_chunks = T[fetch(t) for t in tasks]
#        return 2sum(lnλ_chunks)
#end
