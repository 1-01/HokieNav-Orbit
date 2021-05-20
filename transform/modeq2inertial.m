function [R, varargout] = modeq2inertial(p, f, g, h, k, L, varargin)
% 
% Matt Werner (m.werner@vt.edu) - April 16, 2021
% 
% Calculate the inertial position components using the modified equinoctial
% orbital elements (relative to the ECI coordinate frame) according to
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
%   --- varargin ---
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
%                 R - Inertial position vector consisting of the X, Y, and
%                     Z components relative to the ECI coordinate frame.
%                     Size: N-by-3
%                     Units: km
% 
%   --- varargout ---
% 
%                 V - Inertial velocity vector consisting of the X, Y, and
%                     Z components relative to the ECI coordinate frame.
%                     Size: N-by-3
%                     Units: km/s
% 

% Ensure that the number of inputs are either 6 or 7
narginchk(6, 7)

% Commit some common terms
A2 = h.^2 - k.^2;
s2 = 1 + h.^2 + k.^2;
cosL = cos(L);
sinL = sin(L);
r = p./(1 + f.*cosL + g.*sinL); % Scalar distance

% Compute the inertial position as [X, Y, Z]
R = (r./s2).*[cosL + A2.*cosL + 2*h.*k.*sinL, ...
              sinL - A2.*sinL + 2*h.*k.*cosL, ...
                  2*(h.*sinL - k.*cosL)];

if (nargin == 7)
    % Extract the gravitational parameter as the final input
    GM = varargin{1};
    % Compute the inertial velocity as V = d/dt([X, Y, Z])
    varargout{1} = [+sinL + A2.*sinL - 2*h.*k.*cosL + g - 2*f.*h.*k + A2.*g, ...
                    -cosL + A2.*cosL + 2*h.*k.*sinL - f + 2*g.*h.*k + A2.*f, ...
                            -2*(h.*cosL + k.*sinL + f.*h + g.*k)] ...
                   .*(-sqrt(GM./p)./s2);
end