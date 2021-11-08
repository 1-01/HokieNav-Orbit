function MEE = orbital2modeq(COE)
% 
% Matt Werner (m.werner@vt.edu) - May 25, 2021
% 
% Calculate the modified equinoctial elements (MEE) given the set of
% classical orbital elements (COE).
% 
%    Inputs:
% 
%               COE - The 6 classical (or Keplerian) orbital elements
%                     (COE/KOE) labelled (a, e, I, W, w, f).
%                     Size: 1-by-1 (structure)
%                           6-by-1 (fields)
%                           n-by-1 (each field)
%                     Units: SI (km, radians)
% 
%    Outputs:
% 
%               MEE - The 6 modified equinoctial orbital elements (MEE)
%                     labelled (p, f, g, h, k, L).
%                     Size: 1-by-1 (structure)
%                           6-by-1 (fields)
%                           n-by-1 (each field)
%                     Units: SI (km, radians)
% 

% Transform the classical orbital elements (COE) to modified equinoctial
% elements (MEE)
MEE.p = COE.a.*(1 - COE.e.^2);
MEE.f = COE.e.*cos(COE.w + COE.W);
MEE.g = COE.e.*sin(COE.w + COE.W);
MEE.h = tan(COE.I/2).*cos(COE.W);
MEE.k = tan(COE.I/2).*sin(COE.W);
MEE.L = COE.W + COE.w + COE.f;
