module AuxiliaryFunctions

using HypergeometricFunctions

function w(w_0, z_R, z)
    w = w_0 * ((z^2 + z_R^2)/(z_R^2))^(1/2)
end

function generalized_laguerre_polynomial(n,a,x)
    L_na = binomial(n+a, n)*
        HypergeometricFunctions.M(-n,a+1,x)
end

end
