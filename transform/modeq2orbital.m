function [a, e, I, w, W, v] = modeq2orbital(p, f, g, h, k, L)
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
%  p, f, g, h, k, L - Modified equinoctial orbital elements such that each
%                     element is an N-by-1 vector in the particular order
%                     (p, f, g, h, k, L), where
%                       p - Semiparameter
%                       f - ?
%                       g - ?
%                       h - ?
%                       k - ?
%                       L - True longitude
%                     Size: N-by-1 (vector)
%                     Units: SI
% 
%    Outputs:
% 
%  a, e, I, w, W, v - Classical osculating orbital elements such that each
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
e = sqrt(f.^2 + g.^2);
% Calculate the semimajor axis
a = p./(1 - e.^2);
% Calculate the inclination
I = 2*atan(sqrt(h.^2 + k.^2));
% Calculate the right ascension of the ascending node (RAAN)
W = atan2(k, h);
% Calculate the argument of periapsis
wbar = atan2(g, f);
w = wbar - W;
% Calculate the true anomaly
v = L - wbar;