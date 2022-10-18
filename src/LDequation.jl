## Solve the Lippmann-Schwinger equation in separable potential approximation

function propagator(Λ::Float64,n::Int64)
    return p->p/(p^2+Λ^2)
end

# Create the matrix Iᵢⱼ
function create_I_matrix(terms::Array{Function,1},propagator::Function,Λ::Float64,n::Int64)
    I = zeros(n,n)
    return I
end

## Solve the Lippmann-Schwinger equation 
function LS(model::Model,massgap::Vector{Float64})

end

