clc;
clear;

% v1 = [34355.1604281134;15504.9038569190;-11482.5148748182];
% v2 = [26726.761226;1844.441213;28893.001064];
%

% Rtest = [-0.204845625512718         0.470743465953937        -0.858160159276797
%          0.937666919408811          -0.1570992849457        -0.310000907927605
%         -0.280747249237981         -0.86817072280239        -0.409219474260552];
v1 = [1;1;1];
v2 = [sqrt(3); 0; 0];

a = v1/norm(v1);
b = v2/norm(v2);
%b = [0.743066646034029;0.308764982382103;-0.593730700912727];

crossP = cross(a,b);
dotP = dot(a,b);

v = crossP;
s = norm(v);
c = dotP;

vx = [0, -v(3), v(2);
    v(3), 0, -v(1);
    -v(2), v(1), 0];

R = eye(3) + vx + (vx^2)*1/(1 + c);
%R = 2*(a + b)*(a + b)'/((a + b)'*(a + b)) - eye(3);

% roll = 2.144486425248768;
% pitch = 1.146882858224267;
% yaw = 0.355042854161773;

%R = eul2rotm([yaw, pitch, roll]);

roll = atan2(R(2,1),R(1,1));
%pitch = atan2(-R(3,1),sqrt((R(3,2))^2 + (R(3,3))^2))*180/pi;
pitch = -asin(R(3,1));
yaw = atan2(R(3,2),R(3,3));

disp("roll = " + roll*180/pi + " radians");
disp("pitch = " + pitch*180/pi + " radians");
disp("yaw = " + yaw*180/pi + " radians");

% [R1 R2 R3] = rod2angle(R);
% 
% Rrod = [R1(1), R1(2), R1(3);
%     R2(1), R2(2), R2(3);
%     R3(1), R3(2), R3(3)];
% 
% R = Rrod;
% roll = atan2(R(2,1),R(1,1))*180/pi;
% pitch = atan2(-R(3,1),sqrt((R(3,2))^2 + (R(3,3))^2))*180/pi;
% yaw = atan2(R(3,2),R(3,3))*180/pi;

B = R*a;
A = R'*b;

origin = [0;0;0];

hold on;
plot3([origin(1) a(1)], [origin(2) a(2)], [origin(2) a(3)]);
plot3([origin(1) b(1)], [origin(2) b(2)], [origin(2) b(3)]);
xlabel("x");
ylabel("y");
zlabel("z");
