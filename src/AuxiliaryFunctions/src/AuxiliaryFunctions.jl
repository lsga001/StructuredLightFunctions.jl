module AuxiliaryFunctions

using HypergeometricFunctions

export beam_width, laguerre_polynomial

function beam_width(w_0, z_R, z)
    return w_0 * ((z^2 + z_R^2)/(z_R^2))^(1/2)
end

function laguerre_polynomial(n,a,x)
    return binomial(n+a, n)*
            HypergeometricFunctions.M(-n,a+1,x)
end

end
