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
