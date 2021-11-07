% Initial condition q = q1*i + q2*j + q3*k + q4
IC = [0; 0; 0; 1];

% Runge-Kutta algorithm to propagate the state (quaternion) into the future
% by 100 seconds. The initial condition is specified at time t = 0 and the
% final condition will occur at time t = 100 (seconds)
[t, x] = ode45(@odefunc, [0, 100], IC);

% Calculate the rotation matrix at each instant in time (:,:,t)
R = quat2rotm([x(:,4), x(:,1:3)]);
% Calculate the Euler angles associated with the ZYX Euler sequence
eul = rotm2eul(R(:,:,500))

% Use the 500th rotation matrix as an example
R500 = R(:,:,500);
% Calculate yaw, pitch, and roll
yaw = atan2(R500(2,1), R500(1,1));
pitch = atan2(-R500(3,1), sqrt(R500(3,2)^2 + R500(3,3)^2));
roll = atan2(R500(3,2), R500(3,3));
[yaw, pitch, roll]

for k = 1:size(R,3)
    T = R(:,:,k);
    visualizeRotation(T, [1;1;1])
    pause(0.001)
end

function f = odefunc(t, x)
    % State is x = [q1 q2 q3 q4]
    q1 = x(1);
    q2 = x(2);
    q3 = x(3);
    q4 = x(4);
    
    % Calculate the matrix and angular velocity used within the quaternion
    % dynamics
    Xi = [q4, -q3, q2; q3, q4, -q1; -q2, q1, q4; -q1, -q2, -q3];
    w = angularVelocity(t);
    
    % Form the quaternion dynamics
    f = 0.5*Xi*w;
end

function w = angularVelocity(t)
    % Specify sample angular velocity (in rad/s)
    %
    % In practice, this is a time-varying function determined by
    % instrumentation (the IMU)
    w(1,1) = 0.1;
    w(2,1) = 0.2*sin(t);
    w(3,1) = 0.5*cos(sin(t^2));
end