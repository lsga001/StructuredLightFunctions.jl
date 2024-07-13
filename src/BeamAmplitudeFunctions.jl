module AmplitudeDistributions

include("AuxiliaryFunctions.jl")

""" Laguerre Gaussian Amplitude Distribution """
function LaguerreGaussian(x, y, z, w_0, p, l, n=1, lambda=620e-9)

    k_0 = 2*pi*n/lambda;                # wavenumber [m^-1]
    z_R = pi*n*w_0^2/lambda;            # 
    r, phi = sqrt(x^2 + y^2), atan(y,x) # radius [m], angle[rad]  

    C_lp = sqrt(2*factorial(p)/(pi*factorial(p+abs(l))));

    C_lp*
        (1/(w(w_0, z_R, z)))*
        (r*sqrt(2)/w(w_0, z_R, z))^(abs(l))*
        exp(-r^2 / (w(w_0, z_R, z)^2))*
        generalized_laguerre_polynomial(
                p, abs(l), 2*r^2 / (w(w_0, z_R, z)^2))*
        exp(im*l*phi)*
        exp(im*k_0*r^2*z / (2*(z^2 + z_R^2)))*
        exp(-im*(2*p + abs(l) + 1)*atan(z/z_R))

end

end
