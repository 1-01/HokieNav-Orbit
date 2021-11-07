function M = crossProductMatrix(x)
% 
% Matt Werner (m.werner@vt.edu) - Nov 4, 2021
% 
% Create the cross-product equivalent matrix (M) of a vector (x) such that
% the vector cross-product between two vectors is viewed as a linear
% transformation. That is, given two vectors {u, v}, the cross-product
% equivalent matrix satisfies
% 
%                           cross(u, v) = M*v
% 
% with
%                            _                    _
%                           |   0    -u(3)   u(2)  |
%                           |                      |
%                       M = |  u(3)    0    -u(1)  |.
%                           |                      |
%                           |_ -u(2)  u(1)    0   _|
% 
%    Inputs:
% 
%                 v - Vector whose cross-product equivalent matrix is to be
%                     formed.
%                     Size: 3-by-1 (vector)
%                     Units: -
% 
%    Outputs:
% 
%                 M - The cross-product equivalent matrix of the vector v.
%                     Size: 3-by-3 (matrix)
%                     Units: -
% 

%% Checks
% No checks

%% Calculation
% Form the cross-product equivalent matrix
M = [0, -x(3), x(2); x(3), 0, -x(1); -x(2), x(1), 0];