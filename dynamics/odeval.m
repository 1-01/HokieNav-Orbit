function dxdt = odeval(t, x, pars)
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
%   x = |   p   |,
%       |       |
%       |   f   |
%       |       |
%       |   g   |
%       |       |
%       |   h   |
%       |       |
%       |   k   |
%       |       |
%       |   L   |
%       |__   __|
% 
% where the 6 above quantities are the modified equinoctial orbital
% elements. This set of (non-canonical) parameters eliminate singularities
% encountered by the classical set of orbital parameters (a, e, I, w, W, M)
% when put through the Lagrange Planetary Equations (LPE) (such
% singularities in this case occur when eccentricity e or inclination I
% vanish, indicating that Lagrange's Planetary Equations are invalid for
% circular and/or equatorial orbits). Instead, the modified equinoctial orbital
% elements produce equations of motion that are invalid only for retrograde
% orbits (I = 180 deg), offering a significant improvement upon the LPEs.
% 
% The modified equinoctial orbital elements are indeed related to the
% classical orbital elements via the following transformations.
%             p
%   a = -------------,                    (Semimajor axis)
%        1 - f2 - g2
% 
%         _________
%   e = \/ f2 + g2 ,                      (Orbital eccentricity)
% 
%                _________
%   I = atan(2 \/ h2 + k2 , 1 - h2 - k2), (Orbital inclination),
% 
% 
%   w = atan(gh - fk, fh + gk),           (Arg. of periapsis)
% 
% 
%   W = atan(k, h),                       (RAAN)
% 
% 
%   v = L - atan(g/f),                    (True anomaly)
% 
%  where 
%   - x2 indicates x*x
%   - xy indicates x*y
%   - atan(y,x) is the 4-quadrant inverse tangent
%   - atan(x) is the standard 1-argument inverse tangent.
%  
% The argument of latitude is obtained as
%   u = w + v 
%     = atan(h sin(L) - k cos(L), h cos(L) + k sin(L)).
% 
% These classical orbital parameters are actually the osculating orbital
% elements. Further, the position and velocity vectors expressed in the ECI
% coordinate frame in terms of the modified equinoctial elements are
% 
%            __                                 __
%           |             2                       |
%           |   cos(L) + A  cos(L) + 2hk sin(L)   |
%        r  |                                     |
%   r = --- |             2                       |
%       s*s |   sin(L) - A  sin(L) + 2hk cos(L)   |
%           |                                     |
%           |          2h sin(L) - 2k cos(L)      |
%           |__                                 __|
% 
% 
% and
% 
%            __                                                   __
%           |             2                                   2     |
%           |   sin(L) + A  sin(L) - 2hk cos(L) + g - 2fhk + A  g   |
%       -1  |                                                       | / G M \1/2
%   v = --- |             2                                   2     | | --- |   ,  
%       s*s |  -cos(L) + A  cos(L) + 2hk sin(L) - f + 2ghk + A  f   | \  p  /   
%           |                                                       |
%           |          -2(h cos(L) + k sin(L) + fh + gk)            |
%           |__                                                   __|
% 
% where the constants used in expressing r and v (both vectors in the ECI
% frame) are defined
%         _________
%   A = \/ h2 - k2 ,
%         ____________
%   s = \/1 + h2 + k2 ,
% 
%   r = p/w,
% 
%   w = 1 + f cos(L) + g sin(L).
% 
% Note that the vector r and scalar r are indeed referencing the same
% quantity. Working out the trigonometric simplications, one indeed obtains
% the standard formula for the orbital radius.
%        p      a (1 - e2)
%   r = --- = --------------
%        w     1 + e cos(v)
% 

% Print the integration time
disp("t = " + t)

% Retrieve the standard gravitational parameter of Earth
GM = pars.GM;

% Extract orbital elements for convenience
p = x(1);
f = x(2);
g = x(3);
h = x(4);
k = x(5);
L = x(6);

% Commit some common terms that are used repeatedly in the dynamics
s2 = 1 + h^2 + k^2;
w = 1 + f*cos(L) + g*sin(L);
r = p/w;
we = 7.29211585537707e-05;

% Define the perturbing force F = [R; T; N] in the polar basis situated in
% the orbital plane. This basis has components rhat, thetahat, and zhat,
% where zhat is in the direction of the angular momentum
% 
% Obtain the perturbing gravitational acceleration
[reci, veci] = modeq2inertial(p, f, g, h, k, L, GM);
recf = [cos(we*t), sin(we*t), 0; -sin(we*t), cos(we*t), 0; 0, 0, 1]*reci'; % <-- Incomplete
colatitude = acos(recf(3)/r);
longitude = atan2(recf(2), recf(1));
longitude(longitude < 0) = longitude(longitude < 0) + 2*pi;
g_sphr = [1,  0, 0; 
          0,  0, 1; 
          0, -1, 0] * gravityPerturbation(398600.4418, 6378.1363, ...
                                         pars.degree, pars.order, pars.Cnm, pars.Snm, ...
                                         r, colatitude, longitude);
                                     
% Obtain the perturbing air drag
m = 4; % Spacecraft mass [kg]
[~, rho] = atmosnrlmsise00(1000*norm(recf)-6378136.3, ...
                           pi/2 - colatitude, ...
                           longitude, ...
                           2021, ...
                           168, ...
                           50000+t, ... <-- Sample time of day (UTC)
                           'None', ...  <-- No warning (F107, F107A, APH)
                           'Oxygen');
atmosDensity = rho(6); % Atmospheric density [kg/m3]
S_CD = 0.01*1; % Product of reference area with CD (guessed area and CD) [m2]
vr = sqrt(pars.GM/p)*(f*sin(L) - g*cos(L)); % Radial velocity [km/s]
vt = sqrt(pars.GM/p)*(1 + f*cos(L) + g*sin(L)); % Transverse velocity [km/s]
adrag_sphr = -0.5e3*atmosDensity*S_CD*norm(veci)*[vr; vt; 0]/m; % Aerodynamic drag acceleration [km/s2]

% Total acceleration perturbations
F = g_sphr + adrag_sphr;

% Write the dynamics using modified equinoctial orbital elements.
% These orbital elements produce equations of motion that are singular for
% I = 180 degrees.
dxdt = zeros(6,1);
% The dynamics are Gauss's form of the Lagrange planetary equations in
% modified equinoctial elements obtained from:
%   1. Walker (Celestial Mechanics 36, 1985) pg. 413 Eq. 9
%       - Note the two (2) discrepancies in dg/dt with sources 2. and 3.
%   2. Gondelach & Armellin (arXiv 2018) pg. 40 (Appx A) Eqs. 96-103
%       - "Eqs." 98 and 100 belong to Eqs. 97 and 99
%   3. Read, Younes, Junkins (Tech Science Press 2016) pg. 71 Eqs. 31-36
% and are defined accordingly.
% dxdt(1,1) = dp/dt,             dxdt(4,1) = dh/dt
% dxdt(2,1) = df/dt,             dxdt(5,1) = dk/dt
% dxdt(3,1) = dg/dt,             dxdt(6,1) = dL/dt
dxdt(1,1) = 2*p/w*F(2);
dxdt(2,1) = +F(1)*sin(L) + ((w + 1)*cos(L) + f)*F(2)/w - g*(h*sin(L) - k*cos(L))/w*F(3);
dxdt(3,1) = -F(1)*cos(L) + ((w + 1)*sin(L) + g)*F(2)/w + f*(h*sin(L) - k*cos(L))/w*F(3);
dxdt(4,1) = s2*cos(L)/2/w*F(3);
dxdt(5,1) = s2*sin(L)/2/w*F(3);
dxdt(6,1) = pars.GM*(w/p)^2 + (h*sin(L) - k*cos(L))/w*F(3);
% Include the common factor sqrt(p/GM).
dxdt = dxdt*sqrt(p/pars.GM);