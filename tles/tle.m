function [a, e, I, W, w, f] = tle(file, units)
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
%  a, e, I, W, w, f - The 6 orbital elements.
%                      a: Semimajor axis [km]
%                      e: Eccentricity []
%                      I: Inclination [deg or rad]
%                      W: Right ascension of ascending node (RAAN) [deg or rad]
%                      w: Argument of periapsis [deg or rad]
%                      f: True anomaly [deg or rad]
%                     Size: 1-by-1 (scalar)
%                     Units: See above (for angles, degrees is the default)
%

%% Checks
% No checks

%% Two-line element set
% Open the TLE file and read each of the 3 lines
fid = fopen(file, 'rb');
% L1c = fscanf(fid, '%24c%', 1);
% L2c = fscanf(fid, '%71c%', 1);
% L3c = fscanf(fid, '%71c%', 1);
L1c = fgetl(fid);
L2c = fgetl(fid);
L3c = fgetl(fid);
disp(L1c);
disp(L2c);
disp(L3c);
fprintf('\n\n')
fclose(fid);

% % Open the TLE file and read TLE elements
% fid = fopen(file, 'rb');
% L1 = fscanf(fid, '%24c%*s', 1);
% L2 = fscanf(fid, '%d%6d%*c%5d%*3c%*2f%f%f%5d%*c%*d%5d%*c%*d%d%5d', [1,9]);
% L3 = fscanf(fid, '%d%6d%f%f%f%f%f%f%f', [1,8]);
% fclose(fid);
% 
% % Extract the TLE elements
% epoch = L2(1,4)*86400;          % Epoch Date and Julian Date Fraction
% Db    = L2(1,5);                % Ballistic Coefficient
% I     = L3(1,3);                % Inclination [deg]
% W     = L3(1,4);                % Right Ascension of the Ascending Node [deg]
% e     = str2num(strcat('0.',L3c(27:33)));            % Eccentricity 
% w     = L3(1,6);                % Argument of periapsis [deg]
% M     = L3(1,7);                % Mean anomaly [deg]
% n     = L3(1,8);                % Mean motion [Revs per day]

% Epoch Date and Julian Date Fraction
epoch = str2double(L2c(21:32))*86400;
% Ballistic Coefficient
Db = str2double(strrep(strcat('0.', L2c(54:59), 'e', L2c(60:61)),' ', ''));

% Extract the TLE orbital elements
I = str2double(L3c(9:16));                  % Inclination [deg]
W = str2double(L3c(18:25));                 % Right Ascension of the Ascending Node [deg]
e = str2double(strcat('0.', L3c(27:33)));   % Eccentricity
w = str2double(L3c(35:42));                 % Argument of periapsis [deg]
M = str2double(L3c(44:51));                 % Mean anomaly [deg]
n = str2double(L3c(53:62));                 % Mean motion [revs per day]

%% Orbital elements
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