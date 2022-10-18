## Solving gap equation
function gap(model::Model,initial_guess)
    # extract parameters
    σ, Λ, m0 = model.param.σ, model.param.Λ, model.param.m0

    n = model.n

    # separable approximation of bivariate function
    f = terms(model.pot,Λ,n)

    m(q,ϕ) = m0 + sum(ϕ .* map.(f,q))
    En(q,ϕ) = sqrt(q^2+m(q,ϕ)^2);

    function gap_eqn!(F,ϕ)
        for i=1:length(ϕ)
            F[i] = ϕ[i] - (12V/(π^2))*quadgk(q->q^2*f[i](q)*m(q,ϕ)/En(q,ϕ),0,Inf)[1]
        end
    end
    sol = nlsolve(gap_eqn!, initial_guess, autodiff=:forward)

    return sol.zero
end

export gap