function adrag = dragPerturbation(BStar, density, velocity)
% 
% Matt Werner (m.werner@vt.edu) - June 7 2021
% 
% Calculate the aerodynamic drag experienced by the spacecraft at its
% current position and time. The reference frame of the acceleration is
% identical to that in which expresses the velocity components.
% 
% Note that the result is an acceleration, not a force. This is enacted by
% the BStar (B*) term, which is proportional to 1/m (where m is the
% spacecraft mass).
% 
% Note that Bstar and density contain METERS but the velocity is given in
% terms of KILOMETERS/s. As the distance units differ, an extra factor of
% 1000 is premultiplied when determining the acceleration to ensure its
% units are consistent with the length dimension of velocity (KILOMETERS).
% 
%    Inputs:
% 
%             BStar - (Also B*) Slop term containing information about the
%                     drag coefficient, drag area, and mass. The relation
%                     between these terms and B* is exactly
%                                       B* = CD*A/m,
%                     where CD is the drag coefficient, A is the area, and
%                     m is the mass.
%                     Size: 1-by-1 (scalar)
%                     Units: m2/kg (square meters per kilogram)
% 
%           density - Atmospheric density at the current position and time.
%                     This quantity is determined by an atmosphere model.
%                     Size: 1-by-1 (scalar)
%                     Units: kg/m3 (kilogram per cubic meter)
% 
%          velocity - Spacecraft velocity relative to the atmosphere. The
%                     coordinate frame used to express this quantity may be
%                     any; the resulting acceleration will be expressed in
%                     the same frame.
%                     Size: 3-by-1 (vector)
%                     Units: km/s (kilometers per second)
% 
%    Outputs:
% 
%             adrag - Acceleration experienced by the spacecraft due to
%                     aerodynamic drag. The coordinate frame expressing its
%                     components is identical to the coordinate frame
%                     expressing the components of the given velocity.
%                     Size: 3-by-1 (vector)
%                     Units: km/s2 (kilometers per squared second)
% 

% Calculate the air drag
adrag = -0.5e3*BStar*density*norm(velocity)*velocity;