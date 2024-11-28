@testset "simpson" begin
        # trivial tests
        @test simpson(x -> 1, 1, 2) == 1
        @test simpson(x -> 1, 1, 1) == 0
        # Linear and quadratic test
        # (this shoud be exactly equal)
        @test simpson(x -> x, 0, 1) == 1 / 2
        @test simpson(x -> x^2, 0, 1) == 1 / 3
        # Other integrals test
        # From now on, the error is expected to be
        # of the order h^5 * f⁽⁴⁾(ξ) with h being the domain of
        # integration and ξ some point in the domain. See
        # https://en.wikipedia.org/wiki/Simpson%27s_rule
        @test simpson(x -> x^3, 1, 1.1) ≈ (1.1^4 - 1) / 4 atol = 1e-5
        @test simpson(x -> 1 / x, 0.4, 0.5) ≈ log(5 / 4) atol = 1e-5
        @test simpson(x -> log(x), 1, 1.1) ≈ 1.1log(1.1) - 0.1 atol = 1e-5
end
