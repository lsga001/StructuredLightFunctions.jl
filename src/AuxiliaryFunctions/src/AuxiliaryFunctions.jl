module AuxiliaryFunctions

export cartesian_to_polar, wavenumber,
    rayleigh_length, gaussian_beam_radius, 
    beam_radius_of_curvature,
    laguerre_polynomial

using HypergeometricFunctions


function cartesian_to_polar(x, y)
    return sqrt(x^2 + y^2), atan(y,x)
end

function wavenumber(lambda, n_index=1)
    return 2*pi*n_index/lambda
end

function rayleigh_length(w_0, lambda, n_index)
    return pi*n_index*w_0^2/lambda
end

function gaussian_beam_radius(w_0, z_R, z)
    return w_0 * ((z^2 + z_R^2)/(z_R^2))^(1/2)
end

function beam_radius_of_curvature(z_R, z)
    return z*(1 + (z_R/z)^2)
end

function laguerre_polynomial(n,a,x)
    return binomial(n+a, n)*
            HypergeometricFunctions.M(-n,a+1,x)
end

end
