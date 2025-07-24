# This file contains code for the one-term separable potential following the Mongan paper.

function imagpart_sep_1_T(ksq, potential, mass, λ)
    return -π * mass * λ * sqrt(ksq) * potential(sqrt(ksq))^2 / (2 * (197.32^2))
end

function realpart_sep_1_T(ksq, potential, mass, λ)
    func(q) = q^2 * potential(q)^2 * PrincipalValue(ksq - q^2) / (197.32)^2
    return mass * λ * integrate(func, 0.0, Inf)
end

function phase_shift_sep_1(ksq, potential, mass, λ)
    return 180 * atan(imagpart_sep_1_T(ksq, potential, mass, λ) / (1 - realpart_sep_1_T(ksq, potential, mass, λ))) / π
end
