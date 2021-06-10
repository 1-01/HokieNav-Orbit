function [r, v] = approxEphemeris(t, body)
% 
% Matt Werner (m.werner@vt.edu) - June 9, 2021
% 
% Approximate the "true" ephemeris integrated by JPL (NAIF) by using a best
% fit valid for years 1800 AD to 2050 AD. This routine approximates the
% Keplerian orbital elements of the following bodies:
%   1. Mercury
%   2. Venus
%   3. EM Barycenter
%   4. Mars
%   5. Jupiter
%   6. Saturn
%   7. Uranus
%   8. Neptune
%   9. Pluto
% in their own respective (heliocentric) orbital planes by a simple
% first-order Taylor series. The position and velocity are then
% transformable into the J2000 ecliptic coordinate frame, which is then
% transformable into the J2000 equatorial coordinate frame. 
% 
% The final coordinate frame reporting the position and velocity of the
% body is the J2000 equatorial frame. This frame is approximately
% equivalent to the ICRF. Thus, no transformation is given to the ICRF.
% 
%    Inputs:
% 
%                 t - Julian date of the ephemeris time.
%                     Size: 1-by-1 (scalar)
%                     Units: JD (days)
% 
%              body - Indication as to which planetary body, or bodies, to
%                     evaluate in the approximate ephemeris. The body is
%                     indicated by spelling the name of the planet. Note
%                     that "Earth" will yield the Earth-moon barycenter
%                     position; not exactly Earth itself. Also note,
%                     however, that the barycenter of the Earth-moon system
%                     lies within Earth (~4000km from Earth's center).
%                     Size: 1-by-n (vector)
%                     Units: N/A (string)
% 
%    Outputs:
% 
%                 r - Approximate heliocentric position of the planetary
%                     body.
%                     Size: 3-by-n (scalar)
%                     Units: km (kilometers)
% 
%                 v - Approximate heliocentric velocity of the planetary
%                     body.
%                     Size: 3-by-n (scalar)
%                     Units: km/s (kilometers per second)
% 

% Check that the given ephemeris time is within the allowable time frame
% (1800 - 2050 AD). Here, 
if (t < 2378496.50000)
    error("Ephemeris time may not predate 1800 AD.")
elseif (t > 2469807.50000)
    error("Ephemeris time may not exceed 2050 AD.")
end

% Uppercase the bodies
body = upper(convertCharsToStrings(body));

% Calculate the number of centuries past J2000.0
T = (t - 2451545.0)/36525;

% Approximate the ephemeris
for k = numel(body):-1:1
    % Obtain the Keplerian orbital elements (relative to the mean ecliptic and
    % equinox of J2000.0) for each request body.
    % 
    % a (AU)  - Semimajor axis
    % e (rad) - Eccentricity
    % I (deg) - Inclination
    % L (deg) - Mean longitude
    % B (deg) - Longitude of perihelion
    % W (deg) - Longitude of the ascending node
    % 
    % Note that distance (a) is in terms of astronomical units (AU) and
    % angles are in terms of degrees (deg).
    switch body(k)
        case "MERCURY"
            a = 0.38709927 + 0.00000037*T;
            e = 0.20563593 + 0.00001906*T;
            I = 7.00497902 - 0.00594749*T;
            L = 252.25032350 + 149472.67411175*T;
            B = 77.45779628 + 0.16047689*T;
            W = 48.33076593 - 0.12534081*T;
        case "VENUS"
            a = 0.72333566 + 0.00000390*T;
            e = 0.00677672 - 0.00004107*T;
            I = 3.39467605 - 0.00078890*T;
            L = 181.97909950 + 58517.81538729*T;
            B = 131.60246718 + 0.00268329*T;
            W = 76.67984255 - 0.27769418*T;
        case {"EARTH", "EM", "EARTH-MOON", "EARTH-MOON BARYCENTER"}
            a = 1.00000261 + 0.00000562*T;
            e = 0.01671123 - 0.00004392*T;
            I = -0.00001531 - 0.01294668*T;
            L = 100.46457166 + 35999.37244981*T;
            B = 102.93768193 + 0.32327364*T;
            W = 0;
        case "MARS"
            a = 1.52371034 + 0.00001847*T;
            e = 0.09339410 + 0.00007882*T;
            I = 1.84969142 - 0.00813131*T;
            L = -4.55343205 + 19140.30268499*T;
            B = -23.94362959 + 0.44441088*T;
            W = 49.55953891 - 0.29257343*T;
        case "JUPITER"
            a = 5.20288700 - 0.00011607*T;
            e = 0.04838624 - 0.00013253*T;
            I = 1.30439695 - 0.00183714*T;
            L = 34.39644051 + 3034.74612775*T;
            B = 14.72847983 + 0.21252668*T;
            W = 100.47390909 + 0.20469106*T;
        case "SATURN"
            a = 9.53667594 - 0.00125060*T;
            e = 0.05386179 - 0.00050991*T;
            I = 2.48599187 + 0.00193609*T;
            L = 49.95424423 + 1222.49362201*T;
            B = 92.59887831 - 0.41897216*T;
            W = 113.66242448 - 0.28867794*T;
        case "URANUS"
            a = 19.18916464 - 0.00196176*T;
            e = 0.04725744 - 0.00004397*T;
            I = 0.77263783 - 0.00242939*T;
            L = 313.23810451 + 428.48202785*T;
            B = 170.95427630 + 0.40805281*T;
            W = 74.01692503 + 0.04240589*T;
        case "NEPTUNE"
            a = 30.06992276 + 0.00026291*T;
            e = 0.00859048 + 0.00005105*T;
            I = 1.77004347 + 0.00035372*T;
            L = -55.12002969 + 218.45945325*T;
            B = 44.96476227 - 0.32241464*T;
            W = 131.78422574 - 0.00508664*T;
        case "PLUTO"
            a = 39.48211675 - 0.00031596*T;
            e = 0.24882730 + 0.00005170*T;
            I = 17.14001206 + 0.00004818*T;
            L = 238.92903833 + 145.20780515*T;
            B = 224.06891629 - 0.04062942*T;
            W = 110.30393684 - 0.01183482*T;
    end
    
    % Compute the argument of perihelion, w, and the mean anomaly, M
    w = B - W;
    M = L - B;
    
    % Modulus the mean anomaly to (-180, 180) degrees and solve Kepler's
    % equation for the eccentric anomaly
    M = mod(M, 360);
    M(M > 180) = M(M > 180) - 360;
    E = rad2deg(anomalies(deg2rad(M), e));
    
    % Write the position and velocity relative to the mean ecliptic and
    % equinox of J2000
    rp = a*[cosd(E) - e; sqrt(1 - e^2)*sind(E); 0]; % Orbital plane position
    recl = R3(-W, "degrees")*R1(-I, "degrees")*R3(-w, "degrees")*rp;
    
    % Transform the position and velocity to the J2000 equatorial
    % coordinate frame (the obliquity at J2000 is 23.43928 degrees)
    req = R1(-23.43928, "degrees")*recl;
    r(:,k) = req;
end

% Convert AU to km
r = 149597870.700*r;