function [p, f, g, h, k, L] = orbital2modeq(a, e, I, w, W, v)
% 
% Matt Werner (m.werner@vt.edu) - May 25, 2021
% 
% Calculate the modified equinoctial elements (MEE) given the set of
% classical orbital elements (COE).
% 
%    Inputs:
% 
%                  a - Semimajor axis
%                      Size: n-by-1 (vector)
%                      Units: km (kilometers)
% 
%                  e - Orbital eccentricity
%                      Size: n-by-1 (vector)
%                      Units: - (unitless)
% 
%                  I - Orbital inclination
%                      Size: n-by-1 (vector)
%                      Units: - (radians)
% 
%                  w - Argument of perigee (AoP)
%                      Size: n-by-1 (vector)
%                      Units: - (radians)
% 
%                  W - Right ascension of the ascending node (RAAN)
%                      Size: n-by-1 (vector)
%                      Units: - (radians)
% 
%                  v - True anomaly
%                      Size: n-by-1 (vector)
%                      Units: - (radians)
% 
%    Outputs:
% 
%                  p - Semiparameter
%                      Size: n-by-1 (vector)
%                      Units: km (kilometers)
% 
%                  f - Component of the eccentricity vector
%                      Size: n-by-1 (vector)
%                      Units: - (unitless)
% 
%                  g - Component of the eccentricity vector
%                      Size: n-by-1 (vector)
%                      Units: - (unitless)
% 
%                  h - Orbit orientation parameter
%                      Size: n-by-1 (vector)
%                      Units: - (unitless)
% 
%                  k - Orbit orientation parameter
%                      Size: n-by-1 (vector)
%                      Units: - (unitless)
% 
%                  L - True anomaly
%                      Size: n-by-1 (vector)
%                      Units: - (radians)
% 

% Transform the classical orbital elements to modified equinoctial elements
p = a.*(1 - e.^2);
f = e.*cos(w + W);
g = e.*sin(w + W);
h = tan(I/2).*cos(W);
k = tan(I/2).*sin(W);
L = W + w + v;
