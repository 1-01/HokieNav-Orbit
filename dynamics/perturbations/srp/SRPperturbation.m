function aSRP = SRPperturbation(t, r, A, CR, m)
% 
% Matt Werner (m.werner@vt.edu) - June 9, 2021
% 
% Calculate the perturbative solar radiation pressure experienced by the
% satellite.
% 
% Note that the spacecraft orbit is assumed to be LEO, so Earth's shadow is
% NOT evaluated as a cone, nor are various regions of the shadow considered
% (umbra, penumbra, etc.). The conditions determining if the spacecraft are
% in Earth's shadow is (Wiesel Modern Astrodynamics Ed. 2, Eq. 4.59)
%                    --->   -->
%                     RS  .  r
%       cos(theta) = ---------- > 0     and     r cos(theta) < Req,
%                      RS   r
% where RS is the spacecraft's position relative to the sun, r is the
% spacecraft's position with respect to Earth, and Req is Earth's
% equatorial radius. (Note: Earth is modeled here as a sphere.)
% 
%    Inputs:
% 
%                 t - Julian date of the ephemeris time. This quantity will
%                     be used to evaluate the approximate position of Earth
%                     in its orbit.
%                     Size: 1-by-1 (scalar)
%                     Units: days
% 
%                 r - Earth-centered inertial (ECI) position of the
%                     spacecraft. This quantity will be used to find the
%                     spacecraft's position with respect to the sun.
%                     Size: 3-by-1 (vector)
%                     Units: km (kilometers)
% 
%                 A - Area susceptible to radiation pressure.
%                     Size: 1-by-1 (scalar)
%                     Units: m2 (square meters)
% 
%                CR - Coefficient of reflectivity.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%                 m - Spacecraft mass.
%                     Size: 1-by-1 (scalar)
%                     Units: kg (kilograms)
% 
%    Outputs:
% 
%              aSRP - Acceleration due to solar radiation pressure.
%                     Note: The acceleration is zero if the spacecraft is
%                     in Earth's shadow. As the layers of the shadow are
%                     not considered, there is no smoothing function
%                     provided for the acceleration; the spacecraft is
%                     either in the shadow (0) or not (1).
%                     Size: 3-by-1 (vector)
%                     Units: km/s2 (kilometers per squared seconds)
% 

% Approximate the ephemeris at this time for Earth's position relative to
% the sun. This evaluation returns the Earth's position relative to
% sun-centered inertial coordinates sharing the same orientation as the
% Earth-centered inertial coordinates in which the spacecraft's position is
% expressed.
RE = approxEphemeris(t, 'earth');

% Obtain the spacecraft's position with respect to the sun
RS = RE + r;

% Determine if in Earth's shadow for LEO orbits (simplest case). If so,
% then return no solar radiation pressure
costheta = RS'*r/norm(RS)/norm(r);
if (costheta > 0 && norm(r)*costheta < 6378.1370)
    aSRP = [0;0;0];
    return
end

% Evaluate solar radiation pressure; its units work out to be km/s2.
% - The solar constant (1365 W/m2) is taken from Modern Astrophysics
%   (Carroll, Ostlie) pg. 61.
% - The expression for solar radiation pressure is taken from:
%   1. GMAT 2020a Mathematical Specifications (DRAFT) Eq. 4.77
%       - Note that the units do not result in an acceleration; there is a
%           missing factor of c included in source (2).
%   2. AI-Solutions
%       - https://ai-solutions.com/_help_Files/spacecraft_cr_nanosecond.htm
aSRP = 1e-3*(1365/299792458)*149597870.7^2*CR*A/(m * norm(RS)^3)*RS;