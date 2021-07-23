function HokieNavPath
% 
% Matt Werner (m.werner@vt.edu) - June 6, 2021
% 
% Set up the path temporarily (for this session of Matlab) to simulate the
% orbital dynamics of HokieNav. These paths are added to the top of the
% Matlab path.
% 
% Any modifications to the path made within this file will be automatically 
% reverted upon closing Matlab per Mathworks default behavior.
% 

% Add necessary paths
addpath('check', ...
        genpath('dynamics'), ... Add subfolders
        'math', ...
        'plot/', ...
        'satglobe4e/', ...
        'time', ...
        'tles', ...
        'transform/')
        