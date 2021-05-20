function r_ecf = modeq2ecf(p, f, g, h, k, L)
% 
% Matt Werner (m.werner@vt.edu) - April 16, 2021
% 
% Calculate the geocentric position components using the modified
% equinoctial orbital elements (relative to the ECF coordinate frame)
% according to
%            __                                 __
%           |             2                       |
%           |   cos(L) + A  cos(L) + 2hk sin(L)   |
%        r  |                                     |
%   r = --- |             2                       |.
%       s*s |   sin(L) - A  sin(L) + 2hk cos(L)   |
%           |                                     |
%           |          2h sin(L) - 2k cos(L)      |
%           |__                                 __|
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
%                     Units: SI (km)
% 
%    Outputs:
% 
%             r_ecf - Inertial position vector consisting of the X, Y, and
%                     Z components relative to the ECF coordinate frame.
%                     Size: N-by-3
%                     Units: km
% 

% Ensure that the number of inputs are either 6 or 7
narginchk(6, 7)

% Commit some common terms
A2 = h.^2 - k.^2;
s2 = 1 + h.^2 + k.^2;
cosL = cos(L);
sinL = sin(L);
r = p./(1 + f.*cosL + g.*sinL); % Scalar distance

% Compute the geocentric position as [x, y, z]
r_ecf = (r./s2).*[cosL + A2.*cosL + 2*h.*k.*sinL, ...
                  sinL - A2.*sinL + 2*h.*k.*cosL, ...
                       2*(h.*sinL - k.*cosL)];