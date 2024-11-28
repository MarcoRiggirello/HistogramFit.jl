@testset "chisquare" begin
        f(x, a) = a[1]
        hp = PoissonianBinsModel(
                edges=(1:11,),
                bincounts=[1 for i in 1:10],
                curve=f,
                params_names=(:const,),
                integrator=simpson
        )
        hm = MultinomialBinsModel(
                edges=(1:11,),
                bincounts=[1 for i in 1:10],
                curve=f,
                params_names=(:const,),
                integrator=simpson
        )
        # Poisson test
        @test chisquare(hp, [1]) == 0
        @test chisquare(hp, [2]) ≈ 2 * 10 * (1 - log(2))
        # Multinomial test
        @test chisquare(hm, [1]) == 0
        @test chisquare(hm, [2]) ≈ -20log(2)
end
