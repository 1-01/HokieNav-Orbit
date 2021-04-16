function n = meanMotion(GM, a)
% 
% Matt Werner (m.werner@vt.edu) - April 11, 2021
% 
% Calculate the mean motion. The mean motion is equivalent to the
% expression
%             __
%       n = 2 || / T,
% 
% where T is the orbital period,
%                     3
%       2      __2   a
%      T  =  4 ||  ------
%                    GM
% 
%    Inputs:
% 
%                GM - Gravitational parameter.
%                     Size: 1-by-1 (scalar)
%                     Units: L3/s2 (cubic distance* per squared second)
% 
%                 a - Semimajor axis.
%                     Size: N-by-1 (scalar)
%                     Units: L (distance*)
% 
%    Outputs:
% 
%                 n - Mean motion.
%                     Size: N-by-1 (scalar)
%                     Units: 1/s (distance*)
%                   (L is either m or km)
% 

% Compute the mean motion
n = sqrt(GM./a.^3);