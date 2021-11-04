function visualizeRotation(R, vA)
% 
% Matt Werner (m.werner@vt.edu) - Oct 30, 2021
% 
% Visualize the effect of a rotation matrix R on the vector v whose
% components are given in the inertial frame (vA) to give the same vector v
% with components in the rotated frame (vB = R*vA). Also see how the
% rotated coordinate basis changes.
% 
%    Inputs:
% 
%                 R - Rotation matrix that transforms the components of a
%                     vector from coordinate system A to coordinate system
%                     B. That is, vB = R*vA so that R = R_BA.
%                     Size: 3-by-3 (matrix)
%                     Units: -
% 
%                vA - Vector expressed in the inertial coordinate system A
%                     whose components are to be rotated into coordinate
%                     system B.
%                     Size: 3-by-1 (vector)
%                     Units: -
% 
%    Outputs:
% 
%                   -
% 

%% Plot
% Set scale of figure
scale = 1.2;

% Create figure and define
O = [0;0;0];
iA = [1;0;0];
jA = [0;1;0];
kA = [0;0;1];
% Plot basis in the A frame
plot3([O(1), iA(1)], [O(2), iA(2)], [O(3), iA(3)], 'k-'), hold on, grid on, text(iA(1), iA(2), iA(3), '$x^A$', 'Interpreter', 'latex')
plot3([O(1), jA(1)], [O(2), jA(2)], [O(3), jA(3)], 'k-'), text(jA(1), jA(2), jA(3), '$y^A$', 'Interpreter', 'latex')
plot3([O(1), kA(1)], [O(2), kA(2)], [O(3), kA(3)], 'k-'), text(kA(1), kA(2), kA(3), '$z^A$', 'Interpreter', 'latex')
% Label the axes in frame A as (x, y, z)
xlabel("$x^A$ (black)", 'Interpreter', 'latex')
ylabel("$y^A$ (black)", 'Interpreter', 'latex')
zlabel("$z^A$ (black)", 'Interpreter', 'latex')
% Set the origin to be in the middle
xlim([-scale,scale]), ylim([-scale,scale]), zlim([-scale,scale])
% Visualize the new set of axes (NOTE: This corresponds to R' to rotate the
% basis vectors in frame A)
ipA = R'*iA;
jpA = R'*jA;
kpA = R'*kA;
% Plot the new basis (that was just rotated) in the A frame
plot3([O(1), ipA(1)], [O(2), ipA(2)], [O(3), ipA(3)], 'b-'), text(ipA(1), ipA(2), ipA(3), "$x^B$", 'Interpreter', 'latex')
plot3([O(1), jpA(1)], [O(2), jpA(2)], [O(3), jpA(3)], 'b-'), text(jpA(1), jpA(2), jpA(3), "$y^B$", 'Interpreter', 'latex')
plot3([O(1), kpA(1)], [O(2), kpA(2)], [O(3), kpA(3)], 'b-'), text(kpA(1), kpA(2), kpA(3), "$z^B$", 'Interpreter', 'latex')

% Plot the given v in terms of the A basis
plot3([O(1), vA(1)], [O(2), vA(2)], [O(3), vA(3)], 'r-', 'MarkerSize', 5), text(vA(1), vA(2), vA(3), "$v^A$", 'Interpreter', 'latex')
% Rotate v from frame A to frame B
vB = R*vA;
% Express vB (which is written in terms of frame B components) as a linear
% combination of frame B's axes written with respect to frame A
M = [ipA, jpA, kpA];
% Plot the vector v in the B frame (expressed in A-frame coordinates)
vBinA = M*vB;
plot3([O(1), vBinA(1)], [O(1), vBinA(2)], [O(1), vBinA(3)], 'r-', 'MarkerSize', 5), text(vBinA(1), vBinA(2), vBinA(3), "$v^B$", 'Interpreter', 'latex')

title(sprintf("$v^A = (%1.3f, %1.3f, %1.3f)$ and $v^B = (%1.3f, %1.3f, %1.3f)$", vA(1), vA(2), vA(3), vB(1), vB(2), vB(3)), 'Interpreter', 'latex')
view(112, 30)
hold off