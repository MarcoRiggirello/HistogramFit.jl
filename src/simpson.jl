# Simpson's rule for integral
# Fallback method to perform integration


function simpson(f, a::Number, b::Number)
        h = (b - a) / 2
        f₀ = f(a)
        f₁ = f(a + h)
        f₂ = f(b)
        return h * (f₀ + 4f₁ + f₂) / 3
end


simpson(f, a::AbstractVector{T}, b::AbstractVector{T}) where {T} = simpson(f, a..., b...)


simpson(f, domain) = simpson(f, domain...)
