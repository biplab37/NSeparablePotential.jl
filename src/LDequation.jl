## Solve the Lippmann-Schwinger equation in separable potential approximation

function propagator(Λ::Float64,n::Int64)
    return p->p/(p^2+Λ^2)
end

# Create the matrix Iᵢⱼ
function create_matrix(terms::Array{Function,1},propagator::Function,Λ::Float64,n::Int64)
    J = zeros(n,n)
    return J
end

## Solve the Lippmann-Schwinger equation 
function LS(model::Model,massgap::Vector{Float64})
    return nothing
end

function denom2(J)
    return 1 + J[1, 1]*J[2, 2] - J[1, 1] - J[2, 2] - J[1, 2]*J[2, 1]
end

function denom3(J)
    return J[1, 2]*J[2, 1] + J[1, 3]*J[3, 1] + J[2, 3]*J[3, 2] + J[1, 3]*J[2, 1]*J[3, 2] + J[1, 2]*J[2, 3]*J[3, 1] + J[1, 1]*J[2, 2]*J[3, 3] + J[1, 1] + J[2, 2] + J[3, 3] - 1.0 - J[1, 1]*J[2, 2] - J[1, 1]*J[3, 3] - J[2, 2]*J[3, 3] - J[1, 3]*J[2, 2]*J[3, 1] - J[1, 1]*J[2, 3]*J[3, 2] - J[1, 2]*J[2, 1]*J[3, 3]
end

