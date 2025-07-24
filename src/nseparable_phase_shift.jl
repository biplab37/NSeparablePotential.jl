# This file contains code to calculate the phase shift with n term separable potentials

function propagator_1(k, q)
    return PrincipalValue(k^2 - q^2)
end

function J_matrix_1_real(k::Number, potentials::Vector, signs::Vector, mass)
    dims = length(signs)
    factor = mass / (197.32)^2
    reJ = zeros(dims, dims)

    for i = 1:dims
        for j = 1:dims
            reJ[i, j] = factor * integrate(x -> x^2 * signs[j] * potentials[i](x) * propagator_1(k, x) * potentials[j](x), 0, Inf)
        end
    end
    return reJ
end

function J_matrix_1_imag(k::Number, potentials::Vector, signs::Vector, mass)
    dims = length(signs)
    factor = mass * (-π) / (197.32)^2
    imJ = zeros(dims, dims)

    for i = 1:dims
        for j = 1:dims
            imJ[i, j] = factor * k * signs[j] * potentials[i](k) * potentials[j](k) / 2
        end
    end

    return imJ
end

function phase_shift_T_matrix(k, mass, potentials::Vector, signs::Vector)
    A = J_matrix_1_real(k, potentials, signs, mass)
    B = J_matrix_1_imag(k, potentials, signs, mass)
    inv_IA = inv(I - A)
    C = inv(I - A + B * inv_IA * B)
    re_IJ_inv = C
    im_IJ_inv = -inv_IA * B * C

    repart_mat = I + re_IJ_inv * A - im_IJ_inv * B
    impart_mat = re_IJ_inv * B - im_IJ_inv * A

    gg = map.(potentials, k)
    ggt = transpose(gg .* signs)

    repart = ggt * repart_mat * gg
    impart = ggt * impart_mat * gg

    return 180 * atan(impart / repart) / π
end
