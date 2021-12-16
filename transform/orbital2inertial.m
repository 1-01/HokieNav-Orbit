function ECI = orbital2inertial(COE)
% 
% Matt Werner (m.werner@vt.edu) - April 10, 2021
% 
% Transform 6 osculating orbital elements (a, e, I, w, W, f) to the
% physical, inertial coordinates R = (X, Y, Z) via rotation matrices.
% 
% The required transformation into the inertial frame is (Murray & Dermott)
% 
%  /   \     /                                           \
%  | X |     | cos(W)cos(w + f) - sin(W)sin(w + f)cos(I) |            2
%  |   |     |                                           |    a (1 - e )
%  | Y |  =  | sin(W)cos(w + f) + cos(W)sin(w + f)cos(I) |  --------------,
%  |   |     |                                           |   1 + e cos(f)
%  | Z |     |              sin(w + f)sin(I)             |  
%  \   /     \                                           /  \_____________/
%                                                                  r
% where f is the true anomaly obtained from Kepler's equation.
% 
%    Inputs:
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
%                     Size: 1-by-1 (structure)
%                           6-by-1 (fields)
%                           n-by-1 (each field)
%                     Units: SI (km, radians)
% 
%    Outputs:
% 
%               ECI - Inertial position vector consisting of the X, Y, and
%                     Z components relative to the ECI coordinate frame.
%                     Size: 1-by-1 (structure)
%                           6-by-1 (fields)
%                           n-by-1 (each field)
%                     Units: SI (km, radians)
% 

% Calculate the orbital radius
r = COE.a.*(1 - COE.e.^2)./(1 + COE.e.*cos(COE.f));

% Calculate the inertial position
L = COE.w + COE.f;
ECI.X = (cos(COE.W).*cos(L) - sin(COE.W).*sin(L).*cos(COE.I)).*r;
ECI.Y = (sin(COE.W).*cos(L) + cos(COE.W).*sin(L).*cos(COE.I)).*r;
ECI.Z = sin(L).*sin(COE.I).*r;