function plotTrajectory(x, y, z)
% 
% Matt Werner (m.werner@vt.edu) - June 3, 2021
% 
% Overlay the simulated trajectory overtop of Earth. The Earth model is
% static and does not rotate nor does it have sun lighting; ECF positions
% (x, y, z) are therefore more accurate, though ECI (X, Y, Z) are still
% viable to plot over top of the Earth. In this case, it should simply be
% understood that Earth is shown at only one instance in time (not
% necessarily pertaining to any particular instant in time during the
% trajectory) while the trajectory is shown over a span of time.
% 
%    Inputs:
% 
%           x, y, z - Time-history of the Cartesian position vectors
%                     tracing the spacecraft's trajectory. Note that the
%                     origin of this coordinate system must be
%                     Earth-centered.
%                     Size: n-by-1 (vector)
%                     Units: km (kilometers)
% 
%    Outputs:
% 
%                   -
% 

% Obtain 3D model of Earth
satglobe4e

% Plot trajectory
hold on
plot3(x, y, z, 'Color', '#A2142F')
plot3(x(end), y(end), z(end), 'w.', 'MarkerSize', 5)
hold off