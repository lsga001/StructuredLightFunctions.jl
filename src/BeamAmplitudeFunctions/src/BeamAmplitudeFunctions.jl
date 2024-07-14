module BeamAmplitudeFunctions


##########
module AmplitudeDistributions

export LaguerreGauss, HermiteGauss

using AuxiliaryFunctions

function LaguerreGauss(x, y, z, w_0, p, l, lambda=620e-9, n_index=1)

    k_0 = wavenumber(lambda, n_index);
    z_R = rayleigh_length(w_0, lambda, n_index);
    r, phi = cartesian_to_polar(x, y);

    C_lp = sqrt(2*factorial(p)/(pi*factorial(p+abs(l))));
    R_z = beam_radius_of_curvature(z_R, z);
    w = gaussian_beam_radius(w_0, z_R, z);

    return C_lp*
        (1/w)*
        (r*sqrt(2)/w)^(abs(l))*
        exp(-r^2 / w^2)*
        laguerre_polynomial(p, abs(l), 2*r^2/w^2)*
        exp(im*l*phi)*
        exp(im*k_0*r^2*z / (2*(z^2 + z_R^2)))*
        exp(-im*(2*p + abs(l) + 1)*atan(z/z_R))

end

function HermiteGauss(x, y, z, w_0, m, n, lambda=620e-9, n_index=1)
    k_0 = wavenumber(lambda, n_index);
    z_R = rayleigh_length(w_0, lambda, n_index);
    
    N = m+n;
    C = sqrt(2 / (pi * factorial(n) * factorial(m))) * 2^(-N/2);
    R_z = beam_radius_of_curvature(z_R, z);

    return nothing 
end

end
##########

##########
module OtherFunctions

export test_fun

function test_fun()
    println("This is a test function in the test module!")
end

end
##########

export AmplitudeDistributions

end
