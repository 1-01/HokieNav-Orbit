function [t, COE] = tle(file, units)
% 
% Source:
%   Small Satellites (2021). Two Line Element
%       (https://www.mathworks.com/matlabcentral/fileexchange/39364-two-line-element),
%       MATLAB Central File Exchange. Retrieved July 23, 2021.
% 
% Modified:
%   Matt Werner (m.werner@vt.edu) - July 22, 2021
% 
% Read the NORAD (T)wo-(L)ine (E)lement set (TLE) to extract the 6 orbital
% elements according to the desired angle units (degrees (default) or
% radians).
% 
%    Inputs:
% 
%              file - Name of the file, either absolute or relative to the
%                     main path.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%             units - Specification as to whether angular units should be
%                     given in degrees (default) or radians.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%    Outputs:
% 
%                 t - Time specified as the Julian date and as a datetime
%                     with respect to UTC and the local (ground station)
%                     timezone.
%                     Size: 1-by-1 (structure)
%                           3-by-1 (fields)
%                           1-by-1 (field)
%                     Units: UTC and local
% 
%               TLE - The 6 orbital elements.
%                      a: Semimajor axis [km]
%                      e: Eccentricity []
%                      I: Inclination [deg or rad]
%                      W: Right ascension of ascending node (RAAN) [deg or rad]
%                      w: Argument of periapsis [deg or rad]
%                      f: True anomaly [deg or rad]
%                     Size: 1-by-1 (structure)
%                     Units: See above (for angles, degrees is the default)
%

%% Checks
% No checks

%% Read the two-line element set
% Open the TLE file and read each of the 3 lines
fid = fopen(file, 'rb');
L1c = fgetl(fid);
L2c = fgetl(fid);
L3c = fgetl(fid);
disp(L1c);
disp(L2c);
disp(L3c);
fprintf('\n\n')
fclose(fid);

%% Time
% Get the last 2 digits of the year and adjust it so that it's the stand
% YYYY format according to NORAD standards (00-56 correspond to 2000-2056
% whereas 57-99 correspond to 1957-1999). Assume that the date is in the
% 2000s and adjust it otherwise.
year = str2double(L2c(19:20));
year = 2000 + year;
if (year >= 2057)
    % The TLE is between 1957-1999, so take off 100 years
    year = year - 100;
end
% Number of seconds since Jan 01 of the year
epoch_YTD_s = str2double(L2c(21:32))*86400;

% Convert the year and number of seconds since Jan 01 into a datetime.
% Take off 1 day since Jan 1 is the first day of the year, but exactly 0
% corresponds to Jan 1 00:00:000 UTC.
datestring = datetime(epoch_YTD_s-86400, 'ConvertFrom', 'epochtime', 'Epoch', year+"-01-01");
% Initial time from the TLE
[t.JD, t.UTC, t.local] = defineInitialTime(datestring, 'd-MMMM-yyyy HH:mm:ss.SSS', '-00:00');

%% Orbital elements
% Extract the TLE orbital elements
I = str2double(L3c(9:16));                  % Inclination [deg]
W = str2double(L3c(18:25));                 % Right Ascension of the Ascending Node [deg]
e = str2double(strcat('0.', L3c(27:33)));   % Eccentricity
w = str2double(L3c(35:42));                 % Argument of periapsis [deg]
M = str2double(L3c(44:51));                 % Mean anomaly [deg]
n = str2double(L3c(53:62));                 % Mean motion [revs per day]

% Semi-major axis
a = (398600.4418/(n*2*pi/86400)^2)^(1/3);
% Calculate the true anomaly using Kepler's equation
[~, f] = anomalies(deg2rad(M), e);
f = rad2deg(f);

%% Display
% Print out the 6 orbital elements
fprintf('  a [km]     e    inc [deg]   RAAN [deg]    w [deg]    f [deg] \n ')
fprintf('%4.2f  %4.4f   %4.4f     %4.4f     %4.4f    %4.4f\n\n', [a e I W w f]);

%% Units
% Change the units from degrees to radians if specified
if (nargin == 2)
    switch units
        case {"deg", "degree", "degrees"}
            % Do nothing
        case {"rad", "radian", "radians"}
            % Switch angle units from degrees to radians
            I = deg2rad(I);
            W = deg2rad(W);
            w = deg2rad(w);
            f = deg2rad(f);
        otherwise
            try
                error("Unrecognized unit (" + units + ") for angles (degrees or radians).")
            catch
                error("Provided units must be a string ("""") or character ('') format.")
            end
    end
end

% Assign values
COE.a = a; COE.e = e; COE.I = I;
COE.W = W; COE.w = w; COE.f = f;