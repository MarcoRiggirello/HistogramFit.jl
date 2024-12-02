
function simpson(f, a::Number, b::Number)
        h = (b - a) / 2
        f₀ = f(a)
        f₁ = f(a + h)
        f₂ = f(b)
        return h * (f₀ + 4f₁ + f₂) / 3
end


simpson(f, a::AbstractVector{T}, b::AbstractVector{T}) where {T} = simpson(f, a..., b...)


simpson(f, domain::Tuple) = simpson(f, domain...)

"""
        simpson(f, domain, α)

Implements the [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule)
to compute 1-d integrals numerically.

Its usage is suggested only when the bins are narrow.

!!! todo
        Explain how narrow is defined (compared to what?)

"""
simpson(f, domain, α) = simpson(x -> f(x, α), domain)


"""
        quadgk(f, domain, α)

Computes 1-d integrals using the `GuadGKJL()` method from `Integrals`.
"""
function quadgk(f, domain, α)
        prob = IntegralProblem(f, domain, α)
        return solve(prob, QuadGKJL()).u
end


"""
        hcubature(f, domain, α)

Computes n-d integrals using the `HCubatureJL()` method from `Integrals`.
"""
function hcubature(f, domain, α)
        prob = IntegralProblem(f, domain, α)
        return solve(prob, HCubatureJL()).u
end
