function [a, e, I, w, W, f] = inertial2orbital(rvec, vvec, GM)
% 
% Matt Werner (m.werner@vt.edu) - May 30 2021
% 
% Calculate the Keplerian orbital elements provided the inertial position
% and velocity.
% 
%    Inputs:
% 
%              rvec - Inertial position.
%                     Size: 3-by-1 (vector)
%                     Units: km (kilometers)
% 
%              vvec - Inertial velocity.
%                     Size: 3-by-1 (vector)
%                     Units: km/s (kilometers per second)
% 
%                GM - Central-body gravitational parameter.
%                     Size: 3-by-1 (vector)
%                     Units: km3/s2 (cubic kilometer per squared second)
% 
%    Outputs:
% 
%  a, e, I, w, W, f - The six (6) Keplerian (classical) orbital elements.
%                     Size: 1-by-1 (scalars)
%                     Units: km (kilometers) and - (radians)
% 

% No checks (note that singularities can occur)

% Magnitudes of the provided position and velocity vectors
r = norm(rvec);
v = norm(vvec);

% (Normalized) specific angular momentum
hAMvec = cross(rvec, vvec);
hAMvecUnit = hAMvec/norm(hAMvec);

% Compute the semimajor axis (a)
a = 1/(2/r - v^2/GM);

% Compute the eccentricity (e)
evec = ((vvec'*vvec - GM/sqrt(rvec'*rvec))*rvec - (rvec'*vvec)*vvec)/GM;
e = norm(evec);

% Compute the inclination (I)
I = acos(hAMvecUnit(3));

% Compute the longitude of the ascending node (W)
nvec = cross([0;0;1], hAMvec);
W = acos([1,0,0] * nvec/norm(nvec));
if (nvec(2) < 0), W = 2*pi - W; end

% Compute the argument of perigee (w)
w = acos(evec'*nvec / (e*norm(nvec)));
if (evec(3) < 0), w = 2*pi - w; end

% Compute the true anomaly
f = acos(rvec'*evec / (r*e));
if (rvec'*vvec < 0), f = 2*pi - f; end