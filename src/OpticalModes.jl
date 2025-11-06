module OpticalModes

import FromFile: @from
@from "./utilities/auxiliary_functions.jl" using AuxiliaryFunctions

using SpecialFunctions

export LG, HG, BG

function LG(x, y, w_0, p, l; z=0, lambda=620e-9, n_index=1)
  k_0 = wavenumber(lambda, n_index);
  z_R = rayleigh_length(w_0, lambda, n_index);

  N = abs(l) + 2*p; #combined mode number
  C = sqrt(2*factorial(p)/(pi*factorial(p+abs(l))));
  C_z = beam_curvature(z_R, z);
  w = gaussian_beam_radius(w_0, z_R, z);
  psi = gouy_phase(z_R, z);

  r, phi = cartesian_to_polar(x, y);

  return C*
    (1/w)*
    (r*sqrt(2)/w)^(abs(l))*
    exp(-r^2 / w^2)*
    laguerre_polynomial(p, abs(l), 2*r^2/w^2)*
    exp(im*l*phi)*
    exp(im*k_0*r^2*C_z/2)*
    exp(-im*(N+1)*psi)
end

function HG(x, y, w_0, m, n; z=0, lambda=620e-9, n_index=1)
  #doi: 10.1088/2040-8978/18/5/055001

  k_0 = wavenumber(lambda, n_index);
  z_R = rayleigh_length(w_0, lambda, n_index);

  N = m+n; #combined mode number
  C = 2^(-N/2)*sqrt(2/(pi*factorial(m)*factorial(n)));
  C_z = beam_curvature(z_R, z);
  w = gaussian_beam_radius(w_0, z_R, z);
  psi = gouy_phase(z_R, z);

  return C*
      (1/w)*
      exp(-(x^2 + y^2) / w^2)*
      hermite_polynomial(m, sqrt(2)*x/w)*
      hermite_polynomial(n, sqrt(2)*y/w)*
      exp(-im*k_0*(x^2 + y^2)*C_z / 2)*
      exp(-im*(N+1)*psi)
end

function BHG(x, y, w0, m, n; lambda=620e-9, n_index=1)
  """
  Binary Hermite-Gauss beam.
  """
  return sign(HG(x, y, w0, m, n; z=0, lambda=lambda, n_index=n_index))*
    exp(-(x^2 + y^2)/w0^2)
end

function BG(x, y, w_0, l, kt; z=0, lambda=620e-9, n_index=1)

  k_0 = wavenumber(lambda, n_index);
  z_R = rayleigh_length(w_0, lambda, n_index);

  C_z = beam_curvature(z_R, z);
  w = gaussian_beam_radius(w_0, z_R, z);
  psi = gouy_phase(z_R, z);
 
  r, phi = cartesian_to_polar(x, y);

  return ( w_0/w )*
      exp( im*(k_0 - kt^2/(2*k_0))*z - psi )*
      besselj(l, kt*r/(1+im*z/z_R) )*
      exp( (-1/w^2 + im*k_0*C_z/2) * (r^2 + kt^2*z^2/k_0^2) )*
      exp(im*l*phi)
end

function IG(x, y, w0, q, m, epsilon; z=0, lambda=620e-9, n_index=1)
  k0 = wavenumber(lambda, n_index);
  zR = rayleigh_length(w0, lambda, n_index);

  w = gaussian_beam_radius(w0, zR, z);
  psi = gouy_phase(zR, z);

  zeta, eta = cartesian_to_elliptic(x, y, w*sqrt(epsilon/2));

  return (w0/w)*
    C(i*eta, p,m)
end

end
