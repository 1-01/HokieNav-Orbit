function groundTrack(longitude, latitude, color)
% 
% Matt Werner (m.werner@vt.edu) - June 2 2021
% 
% Plot the ground track of a satellite. The corresponding longitude and
% (geodetic) latitude therefore come from the Earth-fixed position of the
% spacecraft. The longitude is identical to the spherical coordinate
% representation of longitude, but the latitude is geodetic (not
% geocentric).
% 
%    Inputs:
% 
%         longitude - Longitude of the spacecraft, measured positive
%                     counterclockwise from the prime meridian passing
%                     through Greenwich, England.
%                     Size: n-by-1 (vector)
%                     Units: - (rad)
% 
%          latitude - (Geodetic) latitude of the spacecraft, measured
%                     positive towards the North Pole from the equator.
%                     Size: n-by-1 (vector)
%                     Units: - (rad)
% 
%    Outputs:
% 
%                   -
% 

% Determine how many times the longitude crosses over the antimeridian. By
% taking the (central difference) gradient, this effects to determining how
% many times the gradient (of unit spacing) jumps by -180 degrees
gradLongitude = diff(longitude);
jumpIndices = find(abs(gradLongitude) > 6);

% Plot the coastlines of Earth's continents
load coastlines
plot(coastlon, coastlat, 'k'), hold on, grid on, grid minor
xlim([-180,180]), ylim([-90,90]), xticks(-180:30:180), yticks(-90:30:90)
xlabel("Longitude $\lambda$ [deg]", 'interpreter', 'latex')
ylabel("Latitude $\phi $ [deg]", 'interpreter', 'latex')
set(gcf, 'Position', [930, 350, 962, 407])

% Plot the groundtrack of the spacecraft on top of a projection of Earth
if (numel(jumpIndices) > 0)
    for k = 1:1+numel(jumpIndices)
        if (k == 1)
            % First transit. No discontinuity possible yet so plot
            % everything from the first data point to the first crossing of
            % the antimeridian
            transit = 1:jumpIndices(1);
        elseif (k == numel(jumpIndices)+1)
            % Last transit. No discontinuity possible so plot everything
            % from the last crossing to the end
            transit = 1+jumpIndices(k-1):numel(longitude);
        else
            % Intermediate transit. Discontinuities exist on each side of
            % the transit, so plot everything in between
            transit = 1+jumpIndices(k-1):jumpIndices(k);
        end
        % Plot the current continuous transit from
        plot(rad2deg(longitude(transit)), rad2deg(latitude(transit)), 'Color', '#A2142F');
        hold on
    end
else
    % Otherwise, plot the entire transit since a crossing of the
    % antimeridian doesn't occur
    plot(rad2deg(longitude), rad2deg(latitude), 'Color', '#A2142F');
end