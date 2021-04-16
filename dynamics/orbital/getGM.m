function GM = getGM(unit)
% 
% Matt Werner (mwerner@vt.edu) - April 10, 2021
% 
% Retrieve the Earth's gravitational parameter value according to the
% desired SI units.
% 
%    Inputs:
% 
%              unit - Indication for which SI units of length to use in
%                     returning the earth's gravitational parameter.
%                     Size: 1-by-1 (scalar)
%                     Units: - (N/A)
%                       Permissible options are:
%                            1 - Indicates that units are to be returned
%                                such that length is in METERS.
%                         1000 - Indicates that units are to be returned
%                                such that length is in KILOMETERS.
% 
%    Outputs:
% 
%                GM - Earth's gravitational parameter.
%                     Size: 1-by-1 (scalar)
%                     Units: m3/s2  (cubic meters per squared second)
%                             OR
%                            km3/s2 (cubic kilometers per squared second)

% Provide the earth's gravitational constant according to the desired
% indicated units
switch unit
    case 1
        GM = 3.986004418e14;
    case 1000
        GM = 3.986004418e5;
end