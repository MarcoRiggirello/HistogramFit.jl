@testset "chisquare" begin
        f(_, a) = a[1]
        #####################
        # 1D histogram test #
        #####################
        hp = PoissonianBinsModel(
                edges=(1:11,),
                bincounts=ones(10),
                curve=f,
                params_names=(:const,),
                integrator=simpson
        )
        hm = MultinomialBinsModel(
                edges=(1:11,),
                bincounts=ones(10),
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
        #####################
        # nD histogram test #
        #####################
        function cubature(f, d, α)
                prob = IntegralProblem(f, d, α)
                sol = solve(prob, HCubatureJL(); reltol=1e-3, abstol=1e-3)
                return sol.u
        end
        for n in 2:5
                e = Tuple(1:11 for i in 1:n)
                b = ones((10 for i in 1:n)...)
                hp = PoissonianBinsModel(
                        edges=e,
                        bincounts=b,
                        curve=f,
                        params_names=(:const,),
                        integrator=cubature
                )
                hm = MultinomialBinsModel(
                        edges=e,
                        bincounts=b,
                        curve=f,
                        params_names=(:const,),
                        integrator=cubature
                )
                # Poisson test
                @test chisquare(hp, [1]) ≈ 0 atol = 1e-8
                @test chisquare(hp, [2]) ≈ 2 * 10^n * (1 - log(2))
                # Multinomial test
                @test chisquare(hm, [1]) ≈ 0 atol = 1e-8
                @test chisquare(hm, [2]) ≈ -2 * 10^n * log(2)
        end
end
