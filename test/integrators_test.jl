@testset "simpson" begin
        # trivial tests
        @test simpson(x -> 1, 1, 2) == 1
        @test simpson(x -> 1, 1, 1) == 0
        # Linear and quadratic tests
        # (this shoud be exactly equal)
        @test simpson(x -> x, 0, 1) == 1 / 2
        @test simpson(x -> x^2, 0, 1) == 1 / 3
        # Other integrals tests
        # From now on, the error is expected to be
        # of the order h^5 * f⁽⁴⁾(ξ) with h being the domain of
        # integration and ξ some point in the domain. See
        # https://en.wikipedia.org/wiki/Simpson%27s_rule
        @test simpson(x -> x^3, 1, 1.1) ≈ (1.1^4 - 1) / 4 atol = 1e-5
        @test simpson(x -> 1 / x, 0.4, 0.5) ≈ log(5 / 4) atol = 1e-5
        @test simpson(x -> log(x), 1, 1.1) ≈ 1.1log(1.1) - 0.1 atol = 1e-5
end


@testset "quadgk" begin
        α_null = 1
        # trivial tests
        @test HistogramsFit.quadgk((x, p) -> 1, (1, 1), α_null) == 0
        @test HistogramsFit.quadgk((x, p) -> 1, (1, 2), α_null) ≈ 1
        # Linear and quadratic tests
        @test HistogramsFit.quadgk((x, p) -> x, (0, 1), α_null) ≈ 1 / 2
        @test HistogramsFit.quadgk((x, p) -> x^2, (0, 1), α_null) ≈ 1 / 3
        # Other integrals tests
        @test HistogramsFit.quadgk((x, p) -> x^3, (1, 1.1), α_null) ≈ (1.1^4 - 1) / 4
        @test HistogramsFit.quadgk((x, p) -> 1 / x, (0.4, 0.5), α_null) ≈ log(5 / 4)
        @test HistogramsFit.quadgk((x, p) -> log(x), (1, 1.1), α_null) ≈ 1.1log(1.1) - 0.1
        # Wide integrals tests
        @test HistogramsFit.quadgk((x, p) -> sin(x), (0, 8π), α_null) ≈ 0 atol = 1e-8
end


@testset "hcubature" begin
        α_null = 1
        for n in 2:8
                # trivial tests
                domain = ([1 for _ in 1:n], [1 for _ in 1:n])
                @test HistogramsFit.hcubature((x, p) -> 1, domain, α_null) == 0
                domain = ([1 for _ in 1:n], [2 for _ in 1:n])
                @test HistogramsFit.hcubature((x, p) -> 1, domain, α_null) ≈ 1
                # polynomial tests
                domain = ([0 for _ in 1:n], [1 for _ in 1:n])
                for m in 1:4
                        @test HistogramsFit.hcubature((x, p) -> sum(x .^ p), domain, m) ≈ n / (m + 1)
                end
                # antisymmetrical functions tests
                domain = ([-1 for _ in 1:n], [1 for _ in 1:n])
                for m in 1:2:9
                        @test HistogramsFit.hcubature((x, p) -> prod(x .^ p), domain, m) ≈ 0 atol = 1e-8
                        @test HistogramsFit.hcubature((x, p) -> sum(sin.(x) .^ p), domain, m) ≈ 0 atol = 1e-8
                end
        end
end
