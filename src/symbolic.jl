using Symbolics
using LinearAlgebra

n = 4

@variables g[1:n] J[1:n, 1:n] A[1:n]

Id = Matrix{Float64}(I, n, n)

an = inv(Id - J) * g

terms = Symbolics.scalarize(an)

num, den = Symbolics.arguments(Symbolics.value(simplify(terms[end])))

den
