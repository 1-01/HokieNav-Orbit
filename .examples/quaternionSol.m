clear, clc
% Find roots of quaternion polynomials

% Define 3-1-3 Euler sequence rotation (for testing)
R313 = @(th,ph,ps) [cos(th), sin(th), 0; -sin(th), cos(th), 0; 0, 0, 1]*...
                   [1, 0, 0; 0, cos(ph), sin(ph); 0, -sin(ph), cos(ph)]*...
                   [cos(ps), sin(ps), 0; -sin(ps), cos(ps), 0; 0, 0, 1];

% Define the rotation matrix in terms of the quaternion q with...
%   q = q(1) + q(2)i + q(3)j + q(4)k
% (for evaluating the solution)
Rquat = @(q) [1-2*(q(3).^2+q(4).^2), 2*(q(2).*q(3)-q(4).*q(1)), 2*(q(2).*q(4)+q(3).*q(1))
              2*(q(2).*q(3)+q(4).*q(1)), (1-2*(q(2).^2+q(4).^2)), 2*(q(3).*q(4)-q(2).*q(1))
              2*(q(2).*q(4)-q(3).*q(1)), 2*(q(3).*q(4)+q(2).*q(1)), (1-2*(q(2).^2+q(3).^2))]';

% Make up some data. 
% Here, BA represents the unit vector of the B-field as measured from the
% magnetometer and we're supposing that it's equal to the IGRF calculated
% B-field BB (unit normalized) expressed in the IGRF frame. The two are
% therefore related by a rotation matrix which we want to find.
Rtest = R313(2*pi*rand, 2*pi*rand, 2*pi*rand);
BA = rand(3,1);
BB = Rtest*BA;

% Split out components for clarity.
% - frame A is the s/c frame (Bx, By, Bz)
% - frame B is the IGRF frame (Bx', By', Bz')
Bx = BA(1); By = BA(2); Bz = BA(3);
Bxp = BB(1); Byp = BB(2); Bzp = BB(3);

% Write the equations to be solved. They're all quadratic in the quaternion
% q = [qr; qi; qj; qk].
% Note: This function has a Jacobian with determinant 0, so the Jacobian is
%       not invertible
func = @(q) [Rquat(q)*[Bx; By; Bz] - [Bxp; Byp; Bzp]; q'*q - 1];

% Solve for the quaternion
opts = optimoptions('fsolve', 'MaxIterations', 400, 'Display', 'iter');
q = fsolve(func, [1;0;0;0], opts);
% Form the rotation matrix described by this quaternion and calculate the
% B-field of the measurement in the IGRF frame
Rsol = Rquat(q);
BBsol = Rsol*BA;

% Check properties of the results
disp("-------")
disp("Results")
disp("-------")
% props = [norm(BBsol - BB), norm(Rsol*Rsol'-eye(3)), norm(Rsol'*Rsol-eye(3)), det(Rsol)];
disp("|BBsol - BB| = " + norm(BBsol - BB) + " (should be 0)")
disp("|Rsol*Rsol'| = " + norm(Rsol*Rsol'-eye(3)) + " (should be 0)")
disp("|Rsol'*Rsol| = " + norm(Rsol'*Rsol-eye(3)) + " (should be 0)")
disp("det(Rsol) = " + det(Rsol) + " (should be 1)")
disp("|q| = " + norm(q) + " (should be 1)")