@testset "chisquare_statistics" begin
        g(_, a) = a[1]
        f = IntegralFunction(g)
        #####################
        # 1D histogram test #
        #####################
        hp = PoissonianBinsModel(
                edges=(1:11,),
                bincounts=ones(10),
                curve=f,
                params_names=(:const,),
                integrator=QuadGKJL()
        )
        hm = MultinomialBinsModel(
                edges=(1:11,),
                bincounts=ones(10),
                curve=f,
                params_names=(:const,),
                integrator=QuadGKJL()
        )
        # Poisson test
        @test chisquare_statistics([1], hp) == 0
        @test chisquare_statistics([2], hp) ≈ 2 * 10 * (1 - log(2))
        # Multinomial test
        @test chisquare_statistics([1], hm) == 0
        @test chisquare_statistics([2], hm) ≈ -20log(2)
        #####################
        # nD histogram test #
        #####################
        for n in 2:5
                e = Tuple(1:11 for i in 1:n)
                b = ones((10 for i in 1:n)...)
                hp = PoissonianBinsModel(
                        edges=e,
                        bincounts=b,
                        curve=f,
                        params_names=(:const,),
                        integrator=HCubatureJL()
                )
                hm = MultinomialBinsModel(
                        edges=e,
                        bincounts=b,
                        curve=f,
                        params_names=(:const,),
                        integrator=HCubatureJL()
                )
                # Poisson test
                @test chisquare_statistics([1], hp) ≈ 0 atol = 1e-8
                @test chisquare_statistics([2], hp) ≈ 2 * 10^n * (1 - log(2))
                # Multinomial test
                @test chisquare_statistics([1], hm) ≈ 0 atol = 1e-8
                @test chisquare_statistics([2], hm) ≈ -2 * 10^n * log(2)
        end
end
