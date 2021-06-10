function [dxdt, varargout] = odeval(t, x, earth, IC)
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
%   x = |   X   |,
%       |       |
%       |   Y   |
%       |       |
%       |   Z   |
%       |-------|
%       |   .   |
%       |   X   |
%       |       |
%       |   .   |
%       |   Y   |
%       |       |
%       |   .   |
%       |   Z   |
%       |__   __|
% 
% where the 6 above quantities are the inertial components (Cartesian) of
% position and velocity. This set of elements eliminate singularities
% encountered within the dynamics by other sets of elements (e.g.
% classical/Keplerian) within the Lagrange planetary equations (LPE) (such
% singularities in this case occur when eccentricity e or inclination I
% vanish, indicating that Lagrange's planetary equations are invalid for
% circular and/or equatorial orbits). Instead, singularities are
% encountered upon conversion to the orbital elements after the numerical 
% solution has been produced (such singularities can be handled).
% 
% The Cartesian elements include the main GM/r central potential term;
% subsequent perturbations are therefore numerically overwhelmed, though
% their effects are still effected easily with modern computers
% (Schroeder). Because of this choice of state, however, all 6 elements
% vary rapidly (unlike other choices of elements like Keperlian,
% equinoctial, Delaunay, etc.), resulting in suboptimal computational
% efficiency for determining the trajectory.
% 

% Print the integration time
disp("t = " + t)

%% Setup
% Retrieve the standard gravitational parameter of Earth
GM = earth.GM;

% Extract orbital elements for convenience
reci = x(1:3);
veci = x(4:6);

r = sqrt(reci'*reci);
v = sqrt(veci'*veci);

wEarth = [0; 0; 7.2921150e-5];

%% Time and Position
% Calculate the current time and extract the day number (starting from Jan.
% 1 of the current year) and the seconds passed in the day since midnight
% UTC.
tUTC = IC.t.UTC + t/86400;
DOY = day(tUTC, 'dayofyear');
SOD = hour(tUTC)*3600 + minute(tUTC)*60 + second(tUTC);

% Calculate the current Julian date (UTC)
JDUTC = IC.t.JD + t/86400;

% ================== GREENWICH MEAN SIDEREAL TIME (GMST) ==================
% Calculate Greenwich Mean Sidereal Time (GMST)
[~, GMST] = earthRotationAngles(JDUTC);
% Associated rotation matrices (ECI <--> ECF)
Recf_eci = R3(GMST);
Reci_ecf = Recf_eci';

% Calculate the current position in the Earth-centered fixed (ECF)
% coordinate frames. The ECF position is further expressed in standard
% spherical (SPH) coordinates and geodetic/ellipsoidal (GPS) coordinates.
% 
% ============================ ECF POSITION ===============================
% Obtain the Earth-fixed (noninertial) position and break into components
recf = Recf_eci*reci;
xecf = recf(1); yecf = recf(2); zecf = recf(3);

% ============================ SPH POSITION ===============================
[longitudepipi, geocentricLatitude, ~] = cart2sph(xecf, yecf, zecf);
longitude02pi = longitudepipi;
longitude02pi(longitude02pi < 0) = longitude02pi(longitude02pi < 0) + 2*pi;
geocentricColatitude = pi/2 - geocentricLatitude;
% Associated rotation matrix (ECF <-- SPH)
Recf_sph = sphr2ijk(longitude02pi, geocentricColatitude);

% ============================ GPS POSITION ===============================
% Obtain (ellipsoidal) representation of the spacecraft's current position
% with respect to the ECF frame associated with the earth's WGS84
% ellipsoid.
% 
% The WGS84 ellipsoid is defined by equatorial radius and flattening:
%   Req = 6378.1370 km
%     f = 1/298.257223563
% for which its squared eccentricity is 2*f - f*f. The value below can be
% verified that it is indeed the ellipsoid's eccentricity.
[~, geodeticLatitude, GPSh] = geocentric2geodetic(xecf, yecf, zecf, ...
                                             6378.1370, 0.081819190842621);

% ====================== LOCAL SIDEREAL TIME (LST) ========================
% Calculate the Local Sidereal Time (LST)
LST = GMST + longitudepipi;

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
flags = ones(23,1);
flags(2) = 0; flags(9) = -1;
[~, densities] = atmosnrlmsise00(GPSh*1e3, ... Current GPS height (m)
                           rad2deg(geodeticLatitude), ... latitude (deg)
                           rad2deg(longitude02pi), ... longitude (deg)
                           year(tUTC), ... Year
                           DOY, ... Day of the year
                           SOD, ... Seconds for time of day (UTC)
                           flags, ...
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
% Define the perturbing acceleration ap = [X; Y; Z] in the (inertial) ECI
% basis (with corresponding elements I, J, and K) in which:
%   1. I points towards the vernal equinox on Jan 1, 2000 12:00:000 TT
%        (J200 epoch)
%   2. J completes the right-handed triad of unit vectors (J = K x I)
%   3. K points directly through the North pole at the J2000 epoch
% Note: The alignment of K and k (the 3-axis of the ECF coordinate system)
%       is not exactly perfect in reality. As such, the transformation
%       between the ECI and ECF frames is more complicated than a simple
%       3-rotation, but this rotation provides the majority of information
%       transfer between the frames.

% ============================== GRAVITY ==================================
% Obtain the perturbing gravitational acceleration according to EGM2008.
% Note: The values for Earth's gravitational parameter GM and equatorial
%       radius are specified by EGM2008 for compatibility with Geocentric
%       Coordinate Time (IERS TN No. 36 pg. 79).
gp_sph = gravityPerturbation(398600.4418, 6378.1363, ...
                       earth.degree, earth.order, earth.Cnm, earth.Snm, ...
                       r, geocentricColatitude, longitude02pi);
% Switch around the components from the spherical basis to the ECI basis
gp_eci = Reci_ecf*Recf_sph*gp_sph;

% ========================= AERODYNAMIC DRAG ==============================
% Obtain the perturbing aerodynamic acceleration (drag) according to the
% NRLMMSISE00 atmosphere model at the spacecraft's current position and
% time.
mass = 4; % Spacecraft mass [kg]
S_CD = 0.01*2.4; % Product of reference area with CD (guess) [m2]
adrag_eci = dragPerturbation(S_CD/mass, ...
                             altDensity, ...
                             veci - cross(wEarth, reci));

% ======================= SOLAR RADIATION PRESSURE ========================
% Obtain the perturbing solar radiation pressure (SRP) using an approximate
% ephemeris to evaluate Earth's position relative to the sun. This
% approximation comes from NAIF's best fit to the ephemeris through the
% years 1800 to 2050.
aSRP_eci = SRPperturbation(JDUTC, reci, 0.01, 1.8, mass);

%% Nonlinear Dynamics
% Write the dynamics using standard Cartesian components of position and
% velocity in an inertial frame. These orbital elements produce equations
% of motion that are free of singularities but may result in singularities
% when converting to other orbital elements (e.g. classical).
dxdt = zeros(6,1);
% The dynamics are obtained from Newton's laws of motion which include
% conservative and nonconservative perturbing accelerations to the general
% 2-body problem. Expressions for perturbations are obtained from:
%   1. Wiesel (Modern Astrodynamics, 2010) Ch. 4 & 5
%   2. NASA GSFC (GMAT Mathematical Specifications, 2020) Eq. 4.3
% and are defined accordingly.
% dxdt(1,1) = dX/dt,             dxdt(4,1) = dVX/dt
% dxdt(2,1) = dY/dt,             dxdt(5,1) = dVY/dt
% dxdt(3,1) = dZ/dt,             dxdt(6,1) = dVZ/dt
dxdt(1:3) = veci;
dxdt(4:6) = (-GM/r^2)*(reci/r) + gp_eci ...
                               + adrag_eci ...
                               + aSRP_eci;
%% Variable Output
% ============================= STATE DYNAMICS =============================
varargout{1} = dxdt(1);
varargout{2} = dxdt(2);
varargout{3} = dxdt(3);
varargout{4} = dxdt(4);
varargout{5} = dxdt(5);
varargout{6} = dxdt(6);
% ======================== PERTURBING ACCELERATION ========================
varargout{7} = gp_eci';
varargout{8} = adrag_eci';
% ========================== ATMOSPHEREIC DENSITY =========================
varargout{9} = altDensity;
% ================================= TIME ==================================
varargout{10} = GMST;
% ========================== GEODETIC COORDINATES =========================
varargout{11} = longitudepipi;
varargout{12} = geodeticLatitude;
varargout{13} = GPSh;
% ========================== SPHERICAL COORDINATES ========================
varargout{14} = r;
varargout{15} = geocentricColatitude;
varargout{16} = longitude02pi;
% ========================= EARTH-FIXED COORDINATES =======================
varargout{17} = xecf;
varargout{18} = yecf;
varargout{19} = zecf;
% ======================== EARTH-INERTIAL COORDINATES =====================
varargout{20} = x(1);
varargout{21} = x(2);
varargout{22} = x(3);