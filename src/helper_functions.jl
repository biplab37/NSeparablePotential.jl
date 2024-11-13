function find_zero(f::Function, a::Number, b::Number)
    return nlsolve(x -> f(x...), [(a + b) / 2]).zero[1]
end

function integrate(f::Function, a, b)
    return quadgk(f, a, b)[1]
end

function golden_search_maxima(f, a, b, tol=1e-5)
    inv_ϕ = (sqrt(5) - 1) / 2
    while abs(b - a) > tol
        c = b - (b - a) * inv_ϕ
        d = a + (b - a) * inv_ϕ
        if f(c) > f(d)
            b = d
        else
            a = c
        end
    end
    return (b + a) / 2
end

