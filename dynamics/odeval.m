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

%% Setup
% Retrieve the standard gravitational parameter of Earth
GM = pars.GM;

% Extract orbital elements for convenience
p = x(1); f = x(2); g = x(3); h = x(4); k = x(5); L = x(6);

% Commit some common terms that are used repeatedly in the dynamics
s2 = 1 + h^2 + k^2;
w = 1 + f*cos(L) + g*sin(L);
r = p/w;
we = 7.29211585537707e-05;

%% Position
% Calculate the current position in the Earth-centered inertial (ECI) and
% Earth-centered fixed (ECF) coordinate frames. The ECF position is further
% expressed in standard spherical (SPH) coordinates.

% ============================ ECI POSITION ===============================
% Obtain the Earth-centered (inertial) position and break into components
reci = modeq2inertial(p, f, g, h, k, L);
Xeci = reci(1); Yeci = reci(2); Zeci = reci(3);

% ============================ ECF POSITION ===============================
% Obtain the Earth-fixed (noninertial) position and break into components
recf = [cos(we*t), sin(we*t), 0; -sin(we*t), cos(we*t), 0; 0, 0, 1]*reci'; % <-- Incomplete
xecf = recf(1); yecf = recf(2); zecf = recf(3);

% ============================ SPH POSITION ===============================
[longitude, geocentricLatitude, ~] = cart2sph(xecf, yecf, zecf);
longitude(longitude < 0) = longitude(longitude < 0) + 2*pi;
geocentricColatitude = pi/2 - geocentricLatitude;

%% Atmosphere
% Calculate Earth's atmospheric density at altitude for use regarding drag
% using the Naval Research Lab's Mass Spectrometer and Incoherent Scatter
% Radar Exosphere of 2000 (NRLMSISE00) atmospheric model. This model is
% emperically founded through (many) measurements and is valid from Earth's
% surface through the exosphere. It has good history of being used with
% other satellite simulation software and is comparable to other emperical
% models (like JB2008).

% Utilize NRLMSISE00 to calculate atmospheric densities of individual
% species [1/m3] as well as the total mass density [kg/m3] (element 6)
[~, densities] = atmosnrlmsise00(1000*r - 6378136.3, ... <-- Not quite complete
                           geocentricLatitude, ... <-- Requires geodetic latitude
                           longitude, ...
                           2021, ...
                           168, ...
                           50000+t, ... <-- Sample time of day (UTC)
                           'None', ...  <-- No warning (F107, F107A, APH)
                           'Oxygen');
% Extract the atmospheric density [kg/m3]
altDensity = densities(6);

%% Perturbing Accelerations
% Calculate additional accelerations that are not accounted for in the
% standard 2-body solution. These accelerations can include any other
% (physical) effect that is NOT the central Newtonian gravitational
% acceleration.
% 
% Define the perturbing acceleration ap = [R; T; N] in the spherical-like
% basis in which:
%   1. R points radially outwards
%   2. T points transversely to R in the direction of motion (longitudinal)
%   3. N points normal to the orbital plane in the direction of the angular
%        momentum vector.
% The 3-direction is what differentiates this coordinate system from a
% spherical one.

% ============================== GRAVITY ==================================
% Obtain the perturbing gravitational acceleration according to EGM2008.
% Note: The values for Earth's gravitational parameter GM and equatorial
%       radius are specified by EGM2008 for compatibility with Geocentric
%       Coordinate Time (IERS TN No. 36 pg. 79).
g_sphr = gravityPerturbation(398600.4418, 6378.1363, ...
                             pars.degree, pars.order, pars.Cnm, pars.Snm, ...
                             r, geocentricColatitude, longitude);
% Switch around the components from the spherical basis to the RTN basis
g_sphr = [g_sphr(1); g_sphr(3); -g_sphr(2)];

% ========================= AERODYNAMIC DRAG ==============================
mass = 4; % Spacecraft mass [kg]
S_CD = 0.01*1; % Product of reference area with CD (guessed area and CD) [m2]
vr = sqrt(GM/p)*(f*sin(L) - g*cos(L)); % Radial velocity [km/s]
vt = sqrt(GM/p)*(1 + f*cos(L) + g*sin(L)); % Transverse velocity [km/s]
% Obtain the perturbing air drag [km/s2]
adrag_sphr = -0.5e3*altDensity*S_CD*norm([vr; vt; 0])*[vr; vt; 0]/mass;

% =============================== TOTAL ===================================
% Total acceleration perturbations
ap = g_sphr + adrag_sphr;

%% Nonlinear Dynamics
% Write the dynamics using modified equinoctial orbital elements.
% These orbital elements produce equations of motion that are singular for
% I = 180 degrees.
dxdt = zeros(6,1);
% The dynamics are Gauss's form of the Lagrange planetary equations in
% modified equinoctial elements obtained from:
%   1. Walker (Celestial Mechanics 36, 1985) pg. 413 Eq. 9
%       - 1st and 3rd equations are corrected in its errata
%   2. Gondelach & Armellin (arXiv 2018) pg. 40 (Appx A) Eqs. 96-103
%       - "Eqs." 98 and 100 belong to Eqs. 97 and 99
%   3. Read, Younes, Junkins (Tech Science Press 2016) pg. 71 Eqs. 31-36
% and are defined accordingly.
% dxdt(1,1) = dp/dt,             dxdt(4,1) = dh/dt
% dxdt(2,1) = df/dt,             dxdt(5,1) = dk/dt
% dxdt(3,1) = dg/dt,             dxdt(6,1) = dL/dt
dxdt(1,1) = (2*p/w)*ap(2);
dxdt(2,1) = +sin(L)*ap(1) + ((w + 1)*cos(L) + f)*ap(2)/w - g*(h*sin(L) - k*cos(L))*ap(3)/w;
dxdt(3,1) = -cos(L)*ap(1) + ((w + 1)*sin(L) + g)*ap(2)/w + f*(h*sin(L) - k*cos(L))*ap(3)/w;
dxdt(4,1) = s2/(2*w)*cos(L)*ap(3);
dxdt(5,1) = s2/(2*w)*sin(L)*ap(3);
dxdt(6,1) = GM*(w/p)^2 + (h*sin(L) - k*cos(L))/w*ap(3);
% Include the common factor sqrt(p/GM)
dxdt = sqrt(p/GM)*dxdt;