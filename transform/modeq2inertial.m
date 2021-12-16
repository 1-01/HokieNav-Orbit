function ECI = modeq2inertial(MEE, GM)
% 
% Matt Werner (m.werner@vt.edu) - April 16, 2021
% 
% Calculate the inertial position and velocity components using the
% modified equinoctial orbital elements (relative to the ECI coordinate
% frame) according to
%            __                                 __
%           |             2                       |
%           |   cos(L) + A  cos(L) + 2hk sin(L)   |
%        r  |                                     |
%   R = --- |             2                       |.
%       s*s |   sin(L) - A  sin(L) + 2hk cos(L)   |
%           |                                     |
%           |          2h sin(L) - 2k cos(L)      |
%           |__                                 __|
% 
%    Inputs:
% 
%               MEE - Modified equinoctial orbital elements such that each
%                     element is an N-by-1 vector in the particular order
%                     (p, f, g, h, k, L), where
%                       p - Semiparameter
%                       f - Eccentricity component (1 of 2)
%                       g - Eccentricity component (2 of 2)
%                       h - Orientation parameter  (1 of 2)
%                       k - Orientation parameter  (2 of 2)
%                       L - True longitude
%                     Size: N-by-1 (vector)
%                     Units: SI (km)
% 
%                GM - Gravitational parameter used in determining the
%                     intertial velocity. Note that this quantity must be
%                     provided to calculate the inertial velocity,
%                     otherwise the calculation is skipped and a resulting
%                     error with be thrown if the velocity is requested
%                     regardless.
%                     Size: 1-by-1 (scalar)
%                     Units: km3/s2
% 
%    Outputs:
% 
%               ECI - Inertial position/velocity components consisting of
%                     the elements (X, Y, Z, VX, VY, VZ) relative to the
%                     ECI coordinate frame.
%                     Size: 1-by-1 (structure)
%                           6-by-1 (fields)
%                           n-by-1 (each field)
%                     Units: SI (km, radians)
% 

% Commit some common terms
A2 = MEE.h.^2 - MEE.k.^2;
s2 = 1 + MEE.h.^2 + MEE.k.^2;
cosL = cos(MEE.L);
sinL = sin(MEE.L);
r = MEE.p./(1 + MEE.f.*cosL + MEE.g.*sinL); % Scalar distance

% Compute the inertial position as [X, Y, Z]
ECI.X = (r./s2).*(cosL + A2.*cosL + 2*MEE.h.*MEE.k.*sinL);
ECI.Y = (r./s2).*(sinL - A2.*sinL + 2*MEE.h.*MEE.k.*cosL);
ECI.Z = (r./s2).*2*(MEE.h.*sinL - MEE.k.*cosL);

% Compute the inertial velocity as V = d/dt([X, Y, Z])
ECI.VX = -sqrt(GM./MEE.p).*(+sinL + A2.*sinL - 2*MEE.h.*MEE.k.*cosL + MEE.g - 2*MEE.f.*MEE.h.*MEE.k + A2.*MEE.g)./s2;
ECI.VY = -sqrt(GM./MEE.p).*(-cosL + A2.*cosL + 2*MEE.h.*MEE.k.*sinL - MEE.f + 2*MEE.g.*MEE.h.*MEE.k + A2.*MEE.f)./s2;
ECI.VZ = 2*sqrt(GM./MEE.p).*(MEE.h.*cosL + MEE.k.*sinL + MEE.f.*MEE.h + MEE.g.*MEE.k)./s2;