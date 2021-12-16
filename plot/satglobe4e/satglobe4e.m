function satglobe4e
% SATGLOBE4e - Draw an idealized satellite view of earth
%              scaled to match the WGS84 ellipsoid
%
% This file renders a fully manipulatable satellite view
% of earth at a resolution of four pixels per degree, with added
% international political boundaries and gridlines.
% The imagery was obtained from NASA, then postprocessed.
% The globe was rendered using the Matlab Mapping Toolbox.
%
% The Mapping Toolbox is not needed to use this file; however,
% if you have the toolbox, you will be able to use the plot3m
% command to add your own graphics. If not, you can simply use plot3.
%
% In order to save storage space, this m-file loads image
% data from the file satglobe.mat, and then creates the
% graticule mesh itself. This process allows users who
% do not have the Matlab Mapping Toolbox to render the
% figure, but it does take a few moments to compute the
% mesh. Using this trick, the data storage is reduced
% considerably; however, once the figure is rendered, you
% may wish to save it as a regular Matlab figure file
% to increase speed.
%
% NOTE: Try using set(gcf,'renderer','opengl') for best on-screen
%       performance and set(gcf,'renderer','zbuffer') for
%       best printed performance. (This will vary depeding on
%       your computer system, so try experimenting.)
%       The default renderer is zbuffer in order to obtain
%       decent cross-platform performance. If your image
%       seems sub-optimal, switch to OpenGL.
%
% No warranty; use at your own risk.
%
% ver 1.0 Initial writing, Michael Kleder, 2004
% ver 1.1 Scaled to WGS84 ellipsoid, named "satglobe4e"
%         Michael Kleder, Sep 2005
% 

% Load provided data
load('satglobe4', 'c', 's');

a = 6378.137; % WGS84 earth ellipsoid semimajor axis
b = 6356.75231424518; % semiminor axis from WGS84 flattening coefficient

% Spherical coordinate angles
th = repmat((0:.25:180)'*pi/180, [1,1441]);
ph = repmat((-180:.25:180)*pi/180, [721, 1]);

% Update s
s.children(1).properties.XData = a*sin(th).*cos(ph);
s.children(1).properties.YData = a*sin(th).*sin(ph);
s.children(1).properties.ZData = b*cos(th);
s.children(1).properties.CData = double(c)/255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This code adds (th,ph) gridlines & political borders                %
% for n=6:8                                                             %
%    s.children(n).properties.XData = a*s.children(n).properties.XData; %
%    s.children(n).properties.YData = a*s.children(n).properties.YData; %
%    s.children(n).properties.ZData = b*s.children(n).properties.ZData; %
% end                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make the figure
figure;
struct2handle(s,gcf);
axis equal
axis vis3d
set(gcf, 'color', 'k', 'renderer', 'zbuffer', 'inverthardcopy', 'off');
return