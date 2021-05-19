function a = gravityPerturbationScratch(degree, order, radius, latitude, longitude)
% 
% Matt Werner (m.werner@vt.edu) - May 19, 2021
% 
% Calculate the disturbing gravitational acceleration.
% 
% **This function should NOT be used to calculate the disturbing
% gravitational acceleration within simulation as it loads the geopotential
% coefficients upon each call. It should be used only for checking values
% and/or plotting.
% 
%    Inputs:
% 
%            degree - Maximum degree to use within the geopotential.
%                     expansion
%                     Size: 1-by-1 (scalar)
%                     Units: -
% 
%             order - Maximum order to use within the geopotential.
%                     expansion
%                     Size: 1-by-1 (scalar)
%                     Units: -
% 
%            radius - Orbital radius.
%                     Size: 1-by-1 (scalar)
%                     Units: km (kilometers)
% 
%          latitude - Geocentric latitude (measured from the equator).
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%         longitude - Optional: Longitude measured from Greenwich, England
%                     Size: 1-by-n (vector)
%                     Units: - (radians)
% 
%    Outputs:
% 
%                 a - Gravitation acceleration that is not due to the
%                     Newtonian central force component
%                     Size: 1-by-n (vector)
%                     Units: cm/s2 (centimeters per squared second)
% 

% Check inputs
narginchk(5, 5)

% Load the (normalized) geopotential coefficients
[n, m, ~, Cnm, Snm, ~, ~] ...
    = loadGravitationalCoefficients(degree, order, "tide-free");

% Convert geocentric latitude to the colatitude
colatitude = pi/2 - latitude;

% Calculate the disturbing acceleration
for ii = numel(longitude):-1:1
    a(ii,:) = gravityPerturbation(398600.4418, 6378.1363, ...
                                         n, m, Cnm, Snm, ...
                                         radius, colatitude, longitude(ii))';
end

% Convert distance units from km to cm
a = a*1e5;