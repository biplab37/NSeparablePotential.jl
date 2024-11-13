## Solving gap equation
function gap(model::Model, initial_guess::Vector{Float64}, f::Vector)
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

    sol = nlsolve(gap_eqn!, initial_guess, autodiff=:forward, iterations=1000).zero

    result = ResultGap(massgap=sol, dynamicmass=q -> m(q, sol))

    return result
end

function gap(model::Model, initial_guess)
    f = terms((x, y) -> model.pot(x, y, model.param), model.param.Λ, model.n)
    return gap(model, initial_guess, f)
end

function gap(model::Model)
    f = terms((x, y) -> model.pot(x, y, model.param), model.param.Λ, model.n)
    return gap(model, zeros(model.n), f)
end

function scale_pot(scale::Float64, model::Model, initial_guess::Vector{Float64})
    new_model = Model(model.param, (x, y) -> model.pot(x, y, model.param) * scale)
    return gap(new_model, initial_guess)
end

export gap
