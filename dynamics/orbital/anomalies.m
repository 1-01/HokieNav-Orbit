function [E, f, varargout] = anomalies(M, e)
% 
% Matt Werner (m.werner@vt.edu) - April 10, 2021
% 
% Obtain the true anomaly (f) and eccentric anomaly (E) given the mean
% anomaly (M) and orbital eccentricity (e). Additional computations to
% determine the cosine and sine evaluations of the true and eccentric
% anomalies may be completed if requested according to the number of output
% arguments.
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
%                 f - True anomaly.
%                     Size: N-by-1 (vector)
%                     Units: - (radians)
% 
% ----- varargout -----
% 
%              cosE - Cosine of the eccentric anomaly.
%                     Size: N-by-1 (vector)
%                     Units: - (radians)
% 
%              cosf - Cosine of the true anomaly.
%                     Size: N-by-1 (vector)
%                     Units: - (radians)
% 
%              sinf - Sine of the true anomaly.
%                     Size: N-by-1 (vector)
%                     Units: - (radians)
% 

% Obtain the eccentric anomaly from Kepler's equation
E = kepler(M, e);

% Compute the true anomaly and correct it such that the branch cut of
% tangent lies on the positive real axis
f = 2*atan(sqrt((1 + e)./(1 - e)).*tan(E/2));
f(f < 0) = f(f < 0) + 2*pi;

% Perform additional calculations when desired

if (nargout > 2) % At least 3
    % Calculate the cosine of the eccentric anomaly directly
    varargout{1} = cos(E);
    if (nargout > 3) % At least 4
        % Use the previously evaluated cosine of the eccentric anomaly to
        % compute the compute the cosine of the true anomaly indirectly
        varargout{2} = (varargout{1} - e) ./ (1 - e.*varargout{1});
        if (nargout == 5) % At most 5
            % Evaluate the sine of the true anomaly directly
            varargout{3} = sin(f);
        else
            error("Too many outputs.")
        end
    end
end