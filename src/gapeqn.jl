## Solving gap equation
function gap(model::Model, initial_guess::Vector, terms::Vector{Function})
    if length(initial_guess) != model.n
        error("Initial guess must be of length model.n = $(model.n)")
    end
    # extract parameters
    Λ, m0 = model.param.Λ, model.param.m0

    n = model.n

    m(q, ϕ) = m0 + sum(ϕ .* map.(f, q))
    En(q, ϕ) = sqrt(q^2 + m(q, ϕ)^2)

    function gap_eqn!(F, ϕ)
        for i in eachindex(ϕ)
            F[i] = ϕ[i] - (12 / (π^2)) * quadgk(q -> q^2 * f[i](q) * m(q, ϕ) / En(q, ϕ), 0, Inf)[1]
        end
    end

    sol = nlsolve(gap_eqn!, initial_guess, autodiff=:forward).zero

    result = ResultGap(massgap=sol, dynamicmass=q -> m(q, sol))

    return result
end

function gap(model::Model, initial_guess)
    f = terms(model.pot, Λ, n)
    return gap(model, initial_guess, f)
end

function gap(model::Model)
    f = terms(model.pot, Λ, n)
    return gap(model, zeros(model.n), f)
end


export gap
