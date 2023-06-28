## Solve the Lippmann-Schwinger equation in separable potential approximation

function propagator_sigma(p::Float64, P0::Float64, m::ResultGap)
    return 4p^2 / (4(p^2 +  m.dynamicmass(p)^2) - P0^2)
end

# Create the matrix Jᵢⱼ= Σₚ gᵢ(p)G(p) gⱼ(p)
function create_matrix(model::Model, m::ResultGap, f::Vector{Function}, propagator::Function, P0::Float64)
    n = model.n
    J = zeros(n, n)
    for i = 1:n
        for j = 1:i
            J[i,j] = integrate(p -> f[i](p) * f[j](p) * propagator(p, P0, m), 0, Λ)
            J[j,i] = J[i,j]
        end
    end
    return J
end

## Solve the Lippmann-Schwinger equation 
function LS(model::Model)
    intial_guess = zeros(model.n)
    massgap = gap(model, intial_guess)
    return LS(model, massgap)
end

function LS(model::Model, massgap::ResultGap)
    f = terms(model.pot, model.Λ, model.n)
    return LS(model, massgap, f)
end

function LS(model::Model, f::Vector{Function})
    intial_guess = zeros(model.n)
    massgap = gap(model, initial_guess, terms)
    return LS(model, massgap, f)
end

function LS(model::Model, massgap::ResultGap, f::Vector{Function})
    function denominator(P0)
        J = create_matrix(model, massgap, f, P0)
        return denom(J)
    end
    return find_zero(denominator, 0, 2)
end

## Expression for the denominator of the T matrix in separable approximation calculated Mathematica (sample code in `Symbolics.jl` is also given)

function denom(J::Matrix{Float64})
    n = size(J, 1)
    if n == 2
        return denom2(J)
    elseif n == 3
        return denom3(J)
    elseif n == 4
        return denom4(J)
    elseif n == 5
        return denom5(J)
    else
        error("Not implimented yet!")
    end
end

function denom(J::Float64)
    return 1 - J
end

function denom2(J)
    return 1 + J[1, 1] * J[2, 2] - J[1, 1] - J[2, 2] - J[1, 2] * J[2, 1]
end

function denom3(J)
    return J[1, 2] * J[2, 1] + J[1, 3] * J[3, 1] + J[2, 3] * J[3, 2] + J[1, 3] * J[2, 1] * J[3, 2] + J[1, 2] * J[2, 3] * J[3, 1] + J[1, 1] * J[2, 2] * J[3, 3] + J[1, 1] + J[2, 2] + J[3, 3] - 1.0 - J[1, 1] * J[2, 2] - J[1, 1] * J[3, 3] - J[2, 2] * J[3, 3] - J[1, 3] * J[2, 2] * J[3, 1] - J[1, 1] * J[2, 3] * J[3, 2] - J[1, 2] * J[2, 1] * J[3, 3]
end

function denom4(J)
    return (-J[1, 4] + J[1, 4] * J[2, 2] - J[1, 2] * J[2, 4] + J[1, 4] * J[2, 3] * J[3, 2] - J[1, 3] * J[2, 4] * J[3, 2] + J[1, 4] * J[3, 3] - J[1, 4] * J[2, 2] * J[3, 3] + J[1, 2] * J[2, 4] * J[3, 3] - J[1, 3] * J[3, 4] + J[1, 3] * J[2, 2] * J[3, 4] - J[1, 2] * J[2, 3] * J[3, 4]) * J[4, 1] - (J[1, 4] * J[2, 1] + J[2, 4] - J[1, 1] * J[2, 4] + J[1, 4] * J[2, 3] * J[3, 1] - J[1, 3] * J[2, 4] * J[3, 1] - J[1, 4] * J[2, 1] * J[3, 3] - J[2, 4] * J[3, 3] + J[1, 1] * J[2, 4] * J[3, 3] + J[1, 3] * J[2, 1] * J[3, 4] + J[2, 3] * J[3, 4] - J[1, 1] * J[2, 3] * J[3, 4]) * J[4, 2] + (-(J[1, 4] * J[3, 1]) + J[1, 4] * J[2, 2] * J[3, 1] - J[1, 2] * J[2, 4] * J[3, 1] - J[1, 4] * J[2, 1] * J[3, 2] - J[2, 4] * J[3, 2] + J[1, 1] * J[2, 4] * J[3, 2] - J[3, 4] + J[1, 1] * J[3, 4] + J[1, 2] * J[2, 1] * J[3, 4] + J[2, 2] * J[3, 4] - J[1, 1] * J[2, 2] * J[3, 4]) * J[4, 3] + (-((J[1, 3] - J[1, 3] * J[2, 2] + J[1, 2] * J[2, 3]) * J[3, 1]) + (-(J[1, 3] * J[2, 1]) - J[2, 3] + J[1, 1] * J[2, 3]) * J[3, 2] + (1 - J[1, 1] - J[1, 2] * J[2, 1] - J[2, 2] + J[1, 1] * J[2, 2]) * (1 - J[3, 3])) * (1 - J[4, 4])
end


function denom5(J)
    return -(((-(J[1, 5] * J[2, 4]) + J[1, 4] * J[2, 5] + J[1, 5] * J[2, 4] * J[3, 3] - J[1, 4] * J[2, 5] * J[3, 3] - J[1, 5] * J[2, 3] * J[3, 4] + J[1, 3] * J[2, 5] * J[3, 4] + J[1, 4] * J[2, 3] * J[3, 5] - J[1, 3] * J[2, 4] * J[3, 5]) * J[4, 2] - (J[1, 5] * J[2, 4] * J[3, 2] - J[1, 4] * J[2, 5] * J[3, 2] + J[1, 5] * J[3, 4] - J[1, 5] * J[2, 2] * J[3, 4] + J[1, 2] * J[2, 5] * J[3, 4] - J[1, 4] * J[3, 5] + J[1, 4] * J[2, 2] * J[3, 5] - J[1, 2] * J[2, 4] * J[3, 5]) * J[4, 3] - (-J[1, 5] + J[1, 5] * J[2, 2] - J[1, 2] * J[2, 5] + J[1, 5] * J[2, 3] * J[3, 2] - J[1, 3] * J[2, 5] * J[3, 2] + J[1, 5] * J[3, 3] - J[1, 5] * J[2, 2] * J[3, 3] + J[1, 2] * J[2, 5] * J[3, 3] - J[1, 3] * J[3, 5] + J[1, 3] * J[2, 2] * J[3, 5] - J[1, 2] * J[2, 3] * J[3, 5]) * (1 - J[4, 4]) - (-J[1, 4] + J[1, 4] * J[2, 2] - J[1, 2] * J[2, 4] + J[1, 4] * J[2, 3] * J[3, 2] - J[1, 3] * J[2, 4] * J[3, 2] + J[1, 4] * J[3, 3] - J[1, 4] * J[2, 2] * J[3, 3] + J[1, 2] * J[2, 4] * J[3, 3] - J[1, 3] * J[3, 4] + J[1, 3] * J[2, 2] * J[3, 4] - J[1, 2] * J[2, 3] * J[3, 4]) * J[4, 5]) * J[5, 1]) + ((-(J[1, 5] * J[2, 4]) + J[1, 4] * J[2, 5] + J[1, 5] * J[2, 4] * J[3, 3] - J[1, 4] * J[2, 5] * J[3, 3] - J[1, 5] * J[2, 3] * J[3, 4] + J[1, 3] * J[2, 5] * J[3, 4] + J[1, 4] * J[2, 3] * J[3, 5] - J[1, 3] * J[2, 4] * J[3, 5]) * J[4, 1] - (J[1, 5] * J[2, 4] * J[3, 1] - J[1, 4] * J[2, 5] * J[3, 1] - J[1, 5] * J[2, 1] * J[3, 4] - J[2, 5] * J[3, 4] + J[1, 1] * J[2, 5] * J[3, 4] + J[1, 4] * J[2, 1] * J[3, 5] + J[2, 4] * J[3, 5] - J[1, 1] * J[2, 4] * J[3, 5]) * J[4, 3] - (J[1, 5] * J[2, 1] + J[2, 5] - J[1, 1] * J[2, 5] + J[1, 5] * J[2, 3] * J[3, 1] - J[1, 3] * J[2, 5] * J[3, 1] - J[1, 5] * J[2, 1] * J[3, 3] - J[2, 5] * J[3, 3] + J[1, 1] * J[2, 5] * J[3, 3] + J[1, 3] * J[2, 1] * J[3, 5] + J[2, 3] * J[3, 5] - J[1, 1] * J[2, 3] * J[3, 5]) * (1 - J[4, 4]) - (J[1, 4] * J[2, 1] + J[2, 4] - J[1, 1] * J[2, 4] + J[1, 4] * J[2, 3] * J[3, 1] - J[1, 3] * J[2, 4] * J[3, 1] - J[1, 4] * J[2, 1] * J[3, 3] - J[2, 4] * J[3, 3] + J[1, 1] * J[2, 4] * J[3, 3] + J[1, 3] * J[2, 1] * J[3, 4] + J[2, 3] * J[3, 4] - J[1, 1] * J[2, 3] * J[3, 4]) * J[4, 5]) * J[5, 2] - ((J[1, 5] * J[2, 4] * J[3, 2] - J[1, 4] * J[2, 5] * J[3, 2] + J[1, 5] * J[3, 4] - J[1, 5] * J[2, 2] * J[3, 4] + J[1, 2] * J[2, 5] * J[3, 4] - J[1, 4] * J[3, 5] + J[1, 4] * J[2, 2] * J[3, 5] - J[1, 2] * J[2, 4] * J[3, 5]) * J[4, 1] - (J[1, 5] * J[2, 4] * J[3, 1] - J[1, 4] * J[2, 5] * J[3, 1] - J[1, 5] * J[2, 1] * J[3, 4] - J[2, 5] * J[3, 4] + J[1, 1] * J[2, 5] * J[3, 4] + J[1, 4] * J[2, 1] * J[3, 5] + J[2, 4] * J[3, 5] - J[1, 1] * J[2, 4] * J[3, 5]) * J[4, 2] - (-(J[1, 5] * J[3, 1]) + J[1, 5] * J[2, 2] * J[3, 1] - J[1, 2] * J[2, 5] * J[3, 1] - J[1, 5] * J[2, 1] * J[3, 2] - J[2, 5] * J[3, 2] + J[1, 1] * J[2, 5] * J[3, 2] - J[3, 5] + J[1, 1] * J[3, 5] + J[1, 2] * J[2, 1] * J[3, 5] + J[2, 2] * J[3, 5] - J[1, 1] * J[2, 2] * J[3, 5]) * (1 - J[4, 4]) - (-(J[1, 4] * J[3, 1]) + J[1, 4] * J[2, 2] * J[3, 1] - J[1, 2] * J[2, 4] * J[3, 1] - J[1, 4] * J[2, 1] * J[3, 2] - J[2, 4] * J[3, 2] + J[1, 1] * J[2, 4] * J[3, 2] - J[3, 4] + J[1, 1] * J[3, 4] + J[1, 2] * J[2, 1] * J[3, 4] + J[2, 2] * J[3, 4] - J[1, 1] * J[2, 2] * J[3, 4]) * J[4, 5]) * J[5, 3] + (-(J[1, 5] * J[4, 1]) + J[1, 5] * J[2, 2] * J[4, 1] - J[1, 2] * J[2, 5] * J[4, 1] + J[1, 5] * J[2, 3] * J[3, 2] * J[4, 1] - J[1, 3] * J[2, 5] * J[3, 2] * J[4, 1] + J[1, 5] * J[3, 3] * J[4, 1] - J[1, 5] * J[2, 2] * J[3, 3] * J[4, 1] + J[1, 2] * J[2, 5] * J[3, 3] * J[4, 1] - J[1, 3] * J[3, 5] * J[4, 1] + J[1, 3] * J[2, 2] * J[3, 5] * J[4, 1] - J[1, 2] * J[2, 3] * J[3, 5] * J[4, 1] - J[1, 5] * J[2, 1] * J[4, 2] - J[2, 5] * J[4, 2] + J[1, 1] * J[2, 5] * J[4, 2] - J[1, 5] * J[2, 3] * J[3, 1] * J[4, 2] + J[1, 3] * J[2, 5] * J[3, 1] * J[4, 2] + J[1, 5] * J[2, 1] * J[3, 3] * J[4, 2] + J[2, 5] * J[3, 3] * J[4, 2] - J[1, 1] * J[2, 5] * J[3, 3] * J[4, 2] - J[1, 3] * J[2, 1] * J[3, 5] * J[4, 2] - J[2, 3] * J[3, 5] * J[4, 2] + J[1, 1] * J[2, 3] * J[3, 5] * J[4, 2] - J[1, 5] * J[3, 1] * J[4, 3] + J[1, 5] * J[2, 2] * J[3, 1] * J[4, 3] - J[1, 2] * J[2, 5] * J[3, 1] * J[4, 3] - J[1, 5] * J[2, 1] * J[3, 2] * J[4, 3] - J[2, 5] * J[3, 2] * J[4, 3] + J[1, 1] * J[2, 5] * J[3, 2] * J[4, 3] - J[3, 5] * J[4, 3] + J[1, 1] * J[3, 5] * J[4, 3] + J[1, 2] * J[2, 1] * J[3, 5] * J[4, 3] + J[2, 2] * J[3, 5] * J[4, 3] - J[1, 1] * J[2, 2] * J[3, 5] * J[4, 3] - J[4, 5] + J[1, 1] * J[4, 5] + J[1, 2] * J[2, 1] * J[4, 5] + J[2, 2] * J[4, 5] - J[1, 1] * J[2, 2] * J[4, 5] + J[1, 3] * J[3, 1] * J[4, 5] - J[1, 3] * J[2, 2] * J[3, 1] * J[4, 5] + J[1, 2] * J[2, 3] * J[3, 1] * J[4, 5] + J[1, 3] * J[2, 1] * J[3, 2] * J[4, 5] + J[2, 3] * J[3, 2] * J[4, 5] - J[1, 1] * J[2, 3] * J[3, 2] * J[4, 5] + J[3, 3] * J[4, 5] - J[1, 1] * J[3, 3] * J[4, 5] - J[1, 2] * J[2, 1] * J[3, 3] * J[4, 5] - J[2, 2] * J[3, 3] * J[4, 5] + J[1, 1] * J[2, 2] * J[3, 3] * J[4, 5]) * J[5, 4] + ((-J[1, 4] + J[1, 4] * J[2, 2] - J[1, 2] * J[2, 4] + J[1, 4] * J[2, 3] * J[3, 2] - J[1, 3] * J[2, 4] * J[3, 2] + J[1, 4] * J[3, 3] - J[1, 4] * J[2, 2] * J[3, 3] + J[1, 2] * J[2, 4] * J[3, 3] - J[1, 3] * J[3, 4] + J[1, 3] * J[2, 2] * J[3, 4] - J[1, 2] * J[2, 3] * J[3, 4]) * J[4, 1] - (J[1, 4] * J[2, 1] + J[2, 4] - J[1, 1] * J[2, 4] + J[1, 4] * J[2, 3] * J[3, 1] - J[1, 3] * J[2, 4] * J[3, 1] - J[1, 4] * J[2, 1] * J[3, 3] - J[2, 4] * J[3, 3] + J[1, 1] * J[2, 4] * J[3, 3] + J[1, 3] * J[2, 1] * J[3, 4] + J[2, 3] * J[3, 4] - J[1, 1] * J[2, 3] * J[3, 4]) * J[4, 2] + (-(J[1, 4] * J[3, 1]) + J[1, 4] * J[2, 2] * J[3, 1] - J[1, 2] * J[2, 4] * J[3, 1] - J[1, 4] * J[2, 1] * J[3, 2] - J[2, 4] * J[3, 2] + J[1, 1] * J[2, 4] * J[3, 2] - J[3, 4] + J[1, 1] * J[3, 4] + J[1, 2] * J[2, 1] * J[3, 4] + J[2, 2] * J[3, 4] - J[1, 1] * J[2, 2] * J[3, 4]) * J[4, 3] + (-((J[1, 3] - J[1, 3] * J[2, 2] + J[1, 2] * J[2, 3]) * J[3, 1]) + (-(J[1, 3] * J[2, 1]) - J[2, 3] + J[1, 1] * J[2, 3]) * J[3, 2] + (1 - J[1, 1] - J[1, 2] * J[2, 1] - J[2, 2] + J[1, 1] * J[2, 2]) * (1 - J[3, 3])) * (1 - J[4, 4])) * (1 - J[5, 5])
end
