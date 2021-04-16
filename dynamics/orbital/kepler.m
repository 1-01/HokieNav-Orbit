function E = kepler(M, e)
% 
% Matt Werner (m.werner@vt.edu) - April 10, 2021
% 
% Numerically solve Kepler's equation using a root-finding method outlined
% by Danby (1988) that provides quartic convergence.
% 
% Both the eccentric anomaly (E) and mean anomaly (M) are assumed to be
% given in radians, so a reasonable tolerance is chosen as a stopping
% criterion of 1e-10, which corresponds to ~5e-9 degrees.
% 
%    Inputs:
% 
%                 M - Mean anomaly.
%                     Size: N-by-1 (vector)
%                     Units: - (radians)
% 
%                 e - Orbital eccentricity.
%                     Size: N-by-1 (vector)
%                     Units: - (unitless)
% 
%    Outputs:
% 
%                 E - Eccentric anomaly.
%                     Size: N-by-1 (vector)
%                     Units: - (radians)
% 

% Check that the problem is physical and within the assumptions of Kepler's
% equation
if (e < 0)
    error("Orbit is unphysical.")
elseif (e > 1)
    error("Orbit is open.")
end

% Provide an initial guess for the solution. The value 0.85 is recommended
% by Danby (1988).
M = mod(M, 2*pi);
Ei = M + sign(sin(M))*0.85.*e;

% Rewrite Kepler's equation in the form of a root-finding method
%    f(E) = E - e sin(E) - M
err = inf;
while (err > 1e-10)
    % Evaluate the function f(E) and its first 3 derivatives
    fi = Ei - e.*sin(Ei) - M; %    f(E) = E - e sin(E) - M
    fip = 1 - e.*cos(Ei);    %   f'(E) = 1 - e cos(E)
    fipp = -fi + Ei - M;     %  f''(E) = e sin(E)
    fippp = 1 - fip;        % f'''(E) = e cos(E)
    
    % Define delta terms
    di1 = -fi./fip;
    di2 = -fi./(fip + 0.5*di1.*fipp);
    di3 = -fi./(fip + 0.5*di1.*fipp + (1/6)*di2.^2.*fippp);
    
    Eip1 = Ei + di3;
    
    % Calculate the relative error
    err = norm(Eip1 - Ei);
    if (Ei ~= 0)
        err = err/norm(Ei);
    end
    
    % Update the last guess
    Ei = Eip1;
end

% Assign the solution meeting the tolerance of the root-finder
E = Eip1;