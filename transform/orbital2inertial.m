function R = orbital2inertial(a, e, I, w, W, f)
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
%  a, e, I, w, W, f - Classical osculating orbital elements such that each
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
%    Outputs:
% 
%                 R - Inertial position vector consisting of the X, Y, and
%                     Z components relative to the ECI coordinate frame.
%                     Size: N-by-3
%                     Units: SI
% 

% Calculate the orbital radius
r = a.*(1 - e.^2)./(1 + e.*cos(f));

% Calculate the inertial position
X = (cos(W).*cos(w + f) - sin(W).*sin(w + f).*cos(I)).*r;
Y = (sin(W).*cos(w + f) + cos(W).*sin(w + f).*cos(I)).*r;
Z = sin(w + f).*sin(I).*r;

% Calculate the interial position
R = [X, Y, Z];