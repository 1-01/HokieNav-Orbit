function COE = modeq2orbital(MEE)
% 
% Matt Werner (m.werner@vt.edu) - April 16, 2021
% 
% Convert the modified equinoctial orbital elements to the classical
% orbital elements according to the well-known transformations
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
%                     Units: SI
% 
%    Outputs:
% 
%               COE - Classical osculating orbital elements such that each
%                     element is an N-by-1 vector in the particular order
%                     (a, e, I, w, W, f), where
%                       a - Semimajor axis
%                       e - Orbital eccentricity
%                       I - Orbital inclination
%                       w - Argument of periapsis
%                       W - Right ascension of ascending node
%                       f - True anomaly
%                     Size: N-by-1 (vector)
%                     Units: SI
% 

% No checks

% Calculate the eccentricity
COE.e = sqrt(MEE.f.^2 + MEE.g.^2);
% Calculate the semimajor axis
COE.a = MEE.p./(1 - COE.e.^2);
% Calculate the inclination
COE.I = 2*atan(sqrt(MEE.h.^2 + MEE.k.^2));
% Calculate the right ascension of the ascending node (RAAN)
COE.W = atan2(MEE.k, MEE.h);
COE.W(COE.W < 0) = COE.W(COE.W < 0) + 2*pi;
% Calculate the argument of periapsis
COE.w = atan2(MEE.g.*MEE.h - MEE.f.*MEE.k, MEE.f.*MEE.h + MEE.g.*MEE.k);
COE.w(COE.w < 0) = COE.w(COE.w < 0) + 2*pi;
% Calculate the true anomaly
COE.f = MEE.L - (COE.W + COE.w);