%clear, clc, clf

% Set scale of figure
scale = 1.2;

%% Rotation matrices
% Define rotation matrices in terms of Euler angles (313), the quaternion,
% and the Euler axis/angle
R313 = @(th,ph,ps) [cosd(th), sind(th), 0; -sind(th), cosd(th), 0; 0, 0, 1]*...
                   [1, 0, 0; 0, cosd(ph), sind(ph); 0, -sind(ph), cosd(ph)]*...
                   [cosd(ps), sind(ps), 0; -sind(ps), cosd(ps), 0; 0, 0, 1];

% Define the rotation matrix in terms of the quaternion q with...
%   q = q(1) + q(2)i + q(3)j + q(4)k
Rquat = @(q) [1-2*(q(3).^2+q(4).^2), 2*(q(2).*q(3)-q(4).*q(1)), 2*(q(2).*q(4)+q(3).*q(1))
              2*(q(2).*q(3)+q(4).*q(1)), (1-2*(q(2).^2+q(4).^2)), 2*(q(3).*q(4)-q(2).*q(1))
              2*(q(2).*q(4)-q(3).*q(1)), 2*(q(3).*q(4)+q(2).*q(1)), (1-2*(q(2).^2+q(3).^2))]';

% Define the rotation matrix in terms of the quaternion q written in terms
% of the Euler axis, e = (e1, e2, e3), and the Euler angle Phi
%   q = [cos(Phi/2); e1*sin(Phi/2); e2*sin(Phi/2); e3*sin(Phi/2)]
Reaea = @(axis, angle) [(1-2*((axis(2)*sind(angle/2)).^2+(axis(3)*sind(angle/2)).^2)), 2*((axis(1)*sind(angle/2)).*(axis(2)*sind(angle/2))-(axis(3)*sind(angle/2)).*cosd(angle/2)), 2*((axis(1)*sind(angle/2)).*(axis(3)*sind(angle/2))+(axis(2)*sind(angle/2)).*cosd(angle/2))
                       2*((axis(1)*sind(angle/2)).*(axis(2)*sind(angle/2))+(axis(3)*sind(angle/2)).*cosd(angle/2)), (1-2*((axis(1)*sind(angle/2)).^2+(axis(3)*sind(angle/2)).^2)), 2*((axis(2)*sind(angle/2)).*(axis(3)*sind(angle/2))-(axis(1)*sind(angle/2)).*cosd(angle/2))
                       2*((axis(1)*sind(angle/2)).*(axis(3)*sind(angle/2))-(axis(2)*sind(angle/2)).*cosd(angle/2)), 2*((axis(2)*sind(angle/2)).*(axis(3)*sind(angle/2))+(axis(1)*sind(angle/2)).*cosd(angle/2)), (1-2*((axis(1)*sind(angle/2)).^2+(axis(2)*sind(angle/2)).^2))]';

%% Choose a rotation matrix
% Choose a rotation matrix to use
% R = R313(0,0,45);
% R = Rquat([cosd(45/2); 0; 0; 1*sind(45/2)]);
% R = Reaea([0;0;1], 45);

Rtest = roty(atand(1/sqrt(2)))*rotz(45)';
Rquat = [0.577350269189626         0.577350269189626         0.577350269189626
        -0.627510829367687         0.766162995599226        -0.138652166231538
        -0.522395277249846        -0.282242680757336         0.804637958007183];
Rrodr = [0.577350269189626         0.577350269189626         0.577350269189626
        -0.577350269189626         0.788675134594813        -0.211324865405187
        -0.577350269189626        -0.211324865405187         0.788675134594813];

% Pick a vector in the A frame
% vA = rand(3,1);
vA = [1;1;1];
% vA = [34355.1604281134;15504.9038569190;-11482.5148748182]/norm([34355.1604281134;15504.9038569190;-11482.5148748182]);
         
%% Plot
% Create figure and define
figure(1)
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
for k = 1:3
    if (k == 1)
        R = Rtest;
    elseif (k == 2)
        R = Rquat;
    else
        R = Rrodr;
    end
    % Visualize the new set of axes NOTE THAT THIS CORRESPONDS TO R' TO ROTATE
    % THE BASIS VECTORS IN FRAME A
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
end
hold off
title(sprintf("$v^A = (%1.3f, %1.3f, %1.3f)$ and $v^B = (%1.3f, %1.3f, %1.3f)$", vA(1), vA(2), vA(3), vB(1), vB(2), vB(3)), 'Interpreter', 'latex')