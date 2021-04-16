function f = odeval(t, x)
% 
% Matt Werner (m.werner@vt.edu) - April 10, 2021
% 
% Evaluate dynamics of the first-order model
%   .
%   x = f(t, x)
% 
% according to the state
%        __   __
%       |       |
%   x = |   a   |,
%       |       |
%       |   e   |
%       |       |
%       |   I   |
%       |       |
%       |   w   |
%       |       |
%       |   W   |
%       |       |
%       |  tau  |
%       |__   __|
% 
% where the 6 above quantities are the classical orbital parameters defined
%   a - Semimajor axis
%   e - Orbital eccentricity
%   I - Orbital inclination
%   w - Argument of periapsis
%   W - Right ascension of the ascending node (RAAN)
% tau - Time at last periapsis pass
% 

% Print the integration time
disp("t = " + t)

% Provide the standard gravitational parameter of Earth (km3/s2)
GM = getGM(1000);

% Extract orbital elements for convenience
a = x(1);
e = x(2);
I = x(3);
w = x(4);
W = x(5);
tau = x(6);

% Obtain the true anomaly (through the eccentric anomaly via Kepler's
% equation)
n = meanMotion(GM, a);
M = n*(t - tau);
[~, f, cosE, cosf, sinf] = anomalies(M, e);

% Precompute some repetitive terms
coswplusf = cos(w + f);
sinwplusf = sin(w + f);

% Define the perturbing force F = [R; T; N] in the polar basis situated in
% the orbital plane. This basis has components rhat, thetahat, and zhat,
% where zhat is in the direction of the angular momentum
% 
% Total acceleration perturbations
F = [0; 0; 0];

% Write the dynamics (Solar System Dynamics - Murray, Dermott pgs. 54-57)
f = zeros(6,1);
% f(1,1) = da/dt,             f(4,1) = dw/dt
% f(2,1) = de/dt,             f(5,1) = dW/dt
% f(3,1) = dI/dt,             f(6,1) = d(tau)/dt
% 
% The dynamics are written out of order, specifically (1, 2, 3, 5, 4, 6)
% for efficiency and are singular for I = 0 and/or e = 0.
f(1,1) = 2*a^1.5/sqrt(GM*(1 - e^2)) * (F(1)*e*sinf + F(2)*(1 + e*cosf));
f(2,1) = sqrt(a*(1 - e^2)/GM) * (F(1)*sinf + F(2)*(cosf + cosE));
f(3,1) = sqrt(a*(1 - e^2)/GM) * coswplusf / (1 + e*cosf) * F(3);
f(5,1) = sqrt(a*(1 - e^2)/GM) * sinwplusf / (1 + e*cosf) / sin(I) * F(3);
f(4,1) = sqrt(a*(1 - e^2)/GM) * (-F(1)*cosf + F(2)*(2 + e*cosf)*sinf/(1 + e*cosf))/e - f(5,1)*cos(I);
f(6,1) = (3*(tau - t)*sqrt(a)/sqrt(GM*(1 - e^2))*e*sinf + a^2*(1 - e^2)/GM*(2/(1 + e*cosf) - cosf/e))*F(1) ...
       + (3*(tau - t)*sqrt(a)/sqrt(GM*(1 - e^2))*(1 + e*cosf) + a^2*(1 - e^2)/GM*(sinf*(2 + e*cosf)/e/(1 + e*cosf)))*F(2);