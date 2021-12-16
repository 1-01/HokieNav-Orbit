% Define a test rotation matrix that we will try to match
Rtest = roty(atand(1/sqrt(2)))*rotz(45)';

% Make up some data in the body frame. Let the 2nd measurement be a
% perturbation of the first measurement. This corresponds to using the same
% instrument (magnetometer) taking measurements at slightly different
% times. The assumption here is that the satellite is not spinning very
% fast, so it's effectively not rotating over a short time period.
vB1 = [1;1;1];
vB2 = vB1 + [rand; 3*rand; 8*rand]*1e-6;
vB = [vB1, vB2];

% Form the test vectors in the rest frame - they have to be related to vB
% in this way since we're trying to find the test rotation matrix
vA1 = Rtest*vB1;
vA2 = Rtest*vB2;
vA = [vA1, vA2];

% Use Davenport's q-method to determine the best rotation matrix that fits
% between these two measurements (vB) and calculations (vA)
[q, R] = davenportq(vA, vB, [0.5, 0.5]);

% See how we did by switching R and Rtest. (Alternatively, you could see by
% providing R-Rtest - a success shows no blue axes and the representation
% of vB is (0,0,0) since providing the null matrix means that R = Rtest)
visualizeRotation(R-Rtest, vB1)

% Verify that R is indeed a rotation matrix
disp("det(R) = " + det(R) + " (should be 1)")
disp("|R'*R - I| = " + norm(R'*R - eye(3)) + " (should be 0)")
disp("|R*R' - I| = " + norm(R*R' - eye(3)) + " (should be 0)")
disp("|R - Rtest| = " + norm(R - Rtest) + " (should be 0 ideally)")