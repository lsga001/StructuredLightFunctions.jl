module AuxiliaryFunctions

using HypergeometricFunctions

export cartesian_to_polar, wavenumber,
rayleigh_length, gouy_phase,
gaussian_beam_radius,
beam_curvature, 
laguerre_polynomial, hermite_polynomial

function cartesian_to_polar(x, y)
    return sqrt(x^2 + y^2), atan(y,x)
end

function wavenumber(lambda, n_index=1)
    return 2*pi*n_index/lambda
end

function gouy_phase(z_R, z)
    return atan(z/z_R)
end

function rayleigh_length(w_0, lambda, n_index)
    return pi*n_index*w_0^2/lambda
end

function gaussian_beam_radius(w_0, z_R, z)
    return w_0 * ((z^2 + z_R^2)/(z_R^2))^(1/2)
end

function beam_curvature(z_R, z)
    return z / (z^2 + z_R^2)
end

function laguerre_polynomial(n,a,x)
    return binomial(n+a, n)*
            HypergeometricFunctions.M(-n,a+1,x)
end

function hermite_polynomial(n,x)
    return (2x)^n*
            pFq((-n/2,(1-n)/2), (), -1/x^2)
end

end
