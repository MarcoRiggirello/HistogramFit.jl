#=
# Histogram fitting basics

We want to fit Normal-distributed data to find the distribution parameters
=#
using StatsBase, HistogramsFit
using Plots
using ForwardDiff, DifferentiationInterface
using Optimization, OptimizationOptimJL
import Random #hide
Random.seed!(1); # hide

data = randn(10_000)
nothing # hide
#======================================
## Data generation and model definition

We use the StatsBase histogram to store the data. Our histogram is made by
bins of with 0.1 in the range -3:3.
=#
h = Histogram(-3:0.1:3)
append!(h, data)
plot(h, label="data")
#=
Then we define the curve we want to fit the histogram with. The function has
to be of the form `f(data, parameters)` in order to work with `HistogramsFit`
=#
function gaus(x, p)
    N, μ, σ = p
    return N / √(2π * σ^2) * exp(-(x[1] - μ)^2 / (2 * σ^2))
end
#=======================
## Fitting the histogram

Now we want to find the values of the parameters of our model that best fit
our data. The number of events is not fixed by our "experiment" hence we should
use a independent Poissonian statistics for the bin population:
=#
pbm = PoissonianBinsModel(h, gaus, (:N, :μ, :σ))
#=
!!! todo
    improve model display
=#
#=
Now we can use the `chisquare_statistics` method to generate the function to be minimized
=#
# The p parameter is requested by Optimization
# but it is not used by us
initial_params = [100.0, 1.0, 10.0] # far from the answer on purpose
#=
And we are ready to use Optimization! We can mix and match the optimization
algorithm and the automatic differentiation (AD) method we like. Let's use the
robust `BFGS()` and `ForwardDiff()`
=#
optf = OptimizationFunction(chisquare_statistics, Optimization.AutoForwardDiff())
prob = OptimizationProblem(optf, initial_params, pbm)
sol = solve(prob, Optim.BFGS())
#=
!!! todo
    add details about the optimizer options
=#
#============================================
## Parameters error and goodness-of-fit check
============================================#
#=
The number of degrees of freedom is
=#
dof = length(HistogramsFit.bincounts(pbm)) - length(sol.u)
#=
while the χ² statistics in the optimal parameters is
=#
χ² = chisquare_statistics(sol.u, pbm)
#=
perhaps, the χ² is in the range
=#
dof - √dof ≤ χ² ≤ dof + √dof
#=
as expected for a good fit.

We can use AD to compute the covariance matrix as well
=#
covm = inv(hessian(x -> chisquare_statistics(x, pbm), AutoForwardDiff(), sol.u))
#=
Hence we have
=#
for (i, n) in enumerate(pbm.params_names)
    println("$n = $(sol.u[i]) ± $(√covm[i,i])")
end
#=
## Plot the fit
=#
plot(h, label="data")
x = range(-3, 3, length=100)
y = 0.1gaus.(x, (sol.u,)) # to plot correctly we have to multiply by the bin width.
plot!(x, y, label="fit")
