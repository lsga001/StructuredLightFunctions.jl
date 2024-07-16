module AmplitudeDistributions

import FromFile: @from
@from "utilities/auxiliary_functions.jl" using AuxiliaryFunctions

export Gaussian, LaguerreGauss, HermiteGauss, 
InceGaussEven, InceGaussOdd

function Gaussian(x, y, z, w_0, lambda=620e-9,  n_index=1)
    
    k_0 = wavenumber(lambda, n_index);
    z_R = rayleigh_length(w_0, lambda, n_index);
    r, phi = cartesian_to_polar(x,y);

    C_z = beam_curvature(z_R, z);
    w = gaussian_beam_radius(w_0, z_R, z);
    psi = gouy_phase(z_R, z);

    return (w_0/w)*
        exp(-r^2 / w^2)*
        exp(im*k_0*r^2*C_z/2)*
        exp(-im*psi)
end

function LaguerreGauss(x, y, z, w_0, p, l, lambda=620e-9, n_index=1)

    k_0 = wavenumber(lambda, n_index);
    z_R = rayleigh_length(w_0, lambda, n_index);

    N = abs(l) + 2*p; #combined mode number
    C_z = beam_curvature(z_R, z);
    w = gaussian_beam_radius(w_0, z_R, z);
    psi = gouy_phase(z_R, z);
    
    r, phi = cartesian_to_polar(x, y);

    C = sqrt(2*factorial(p)/(pi*factorial(p+abs(l)))); # normalization constant
    
    return C*
        (1/w)*
        (r*sqrt(2)/w)^(abs(l))*
        exp(-r^2 / w^2)*
        laguerre_polynomial(p, abs(l), 2*r^2/w^2)*
        exp(im*l*phi)*
        exp(im*k_0*r^2*C_z/2)*
        exp(-im*(N+1)*psi)

end

function HermiteGauss(x, y, z, w_0, m, n, lambda=620e-9, n_index=1)
    #doi: 10.1088/2040-8978/18/5/055001
    
    k_0 = wavenumber(lambda, n_index);
    z_R = rayleigh_length(w_0, lambda, n_index);
    
    N = m+n; #combined mode number
    C_z = beam_curvature(z_R, z);
    w = gaussian_beam_radius(w_0, z_R, z);
    psi = gouy_phase(z_R, z);
    
    C = 2^(-N/2)*sqrt(2/(pi*factorial(m)*factorial(n))); # normalization constant

    return C*
        (1/w)*
        exp(-(x^2 + y^2) / w^2)*
        hermite_polynomial(m, sqrt(2)*x/w)*
        hermite_polynomial(n, sqrt(2)*y/w)*
        exp(-im*k_0*(x^2 + y^2)*C_z / 2)*
        exp(-im*(N+1)*psi)
end

function InceGaussEven(x, y, z, w_0, p, m, epsilon, lambda=620e-9, n_index=1)

    k_0 = wavenumber(lambda, n_index);
    z_R = rayleigh_length(w_0, lambda, n_index);

    C_z = beam_curvature(z_R, z);
    w = gaussian_beam_radius(w_0, z_R, z);
    psi = gouy_phase(z_R, z);

    zeta, eta = cartesian_to_elliptical(w*sqrt(epsilon/2), x, y);

    C = 1; # normalization constant
    
    return C*
        (w_0/w)*
        ince_polynomial_even(p, m, im*zeta, epsilon)*
        ince_polynomial_even(p, m, eta, epsilon)*
        exp(-(x^2+y^2) / w^2)*
        exp(im*k_0*(x^2+y^2)*C_z / 2)*
        exp(-im*(p+1)*psi)
end

function InceGaussOdd(x, y, z, w_0, p, m, epsilon, lambda=620e-9, n_index=1)

    k_0 = wavenumber(lambda, n_index);
    z_R = rayleigh_length(w_0, lambda, n_index);

    C_z = beam_curvature(z_R, z);
    w = gaussian_beam_radius(w_0, z_R, z);
    psi = gouy_phase(z_R, z);

    zeta, eta = cartesian_to_elliptical(w*sqrt(epsilon/2), x, y);
    
    S = 1; # normalization constant

    return S*
        (w_0/w)*
        ince_polynomial_odd(p, m, im*zeta, epsilon)*
        ince_polynomial_odd(p, m, eta, epsilon)*
        exp(-(x^2+y^2) / w^2)*
        exp(im*k_0*(x^2+y^2)*C_z / 2)*
        exp(-im*(p+1)*psi)
end

end
