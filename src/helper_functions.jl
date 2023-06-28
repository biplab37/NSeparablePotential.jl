function find_zero(f::Function, a::Number, b::Number)
    return nlsolve(x->f(x...), [a, b], autodiff=:forward).zero[1]
end

function integrate(f::Function, a, b)
    return quadgk(f, a, b)[1]
end
