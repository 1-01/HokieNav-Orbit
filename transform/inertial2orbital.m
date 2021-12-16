function COE = inertial2orbital(ECI, GM)
% 
% Matt Werner (m.werner@vt.edu) - May 30 2021
% 
% Calculate the Keplerian orbital elements provided the inertial position
% and velocity.
% 
%    Inputs:
% 
%               ECI - Inertial position and velocity, labelled by the
%                     Cartesian coordinates (X, Y, Z, VX, VY, VZ).
%                     Size: 1-by-1 (structure)
%                           6-by-1 (fields)
%                           1-by-1 (each field)
%                     Units: SI (km, radians)
% 
%                GM - Central-body gravitational parameter.
%                     Size: 1-by-1 (scalar)
%                     Units: km3/s2 (cubic kilometer per squared second)
% 
%    Outputs:
% 
%               COE - The six (6) Keplerian (classical) orbital elements
%                     labelled (a, e, I, W, w, f).
%                     Size: 1-by-1 (structure)
%                           6-by-1 (fields)
%                           1-by-1 (each field)
%                     Units: SI (km, radians)
% 

% No checks (note that singularities can occur)

rvec = [ECI.X, ECI.Y, ECI.Z];
vvec = [ECI.VX, ECI.VY, ECI.VZ];

% Magnitudes of the provided position and velocity vectors
r = vecnorm(rvec, 2, 2);
v = vecnorm(vvec, 2, 2);

% (Normalized) specific angular momentum
hAMvec = cross(rvec, vvec);
hAMvecUnit = hAMvec./vecnorm(hAMvec,2,2);

% Compute the semimajor axis (a)
a = 1./(2./r - v.^2/GM);

% Compute the eccentricity (e)
evec = ((v.^2 - GM./r).*rvec - sum(rvec.*vvec,2).*vvec)/GM;
e = vecnorm(evec,2,2);

% Compute the inclination (I)
I = acos(hAMvecUnit(:,3));

% Compute the longitude of the ascending node (W)
khat = [zeros(numel(r), 2), ones(numel(r), 1)];
nvec = cross(khat, hAMvec);
W = acos(nvec(:,1)./vecnorm(nvec,2,2));
for k = 1:numel(W)
    if (nvec(k,2) < 0), W(k) = 2*pi - W(k); end
end

% Compute the argument of perigee (w)
w = acos(sum(evec.*nvec,2) ./ (e.*vecnorm(nvec,2,2)));
for k = 1:numel(w)
    if (evec(k,3) < 0), w(k) = 2*pi - w(k); end
end

% Compute the true anomaly
f = acos(sum(rvec.*evec,2) ./ (r.*e));
for k = 1:numel(f)
    if (sum(rvec(k,:).*vvec(k,:)) < 0), f(k) = 2*pi - f(k); end
end

% Assign
COE.a = a; COE.e = e; COE.I = I;
COE.W = W; COE.w = w; COE.f = f;