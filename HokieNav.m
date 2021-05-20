% HokieNav Orbital Simulation
% Matt Werner (m.werner@vt.edu) - April 8, 2021
clear, clc

%% Inputs
% Gravitational parameter (km3/s2)
gravity.GM = getGM(1000);

[gravity.degree, gravity.order, ...
 ~, gravity.Cnm, gravity.Snm, ~, ~] ...
    = loadGravitationalCoefficients(10, 10, "tide-free");

% Initial conditions read from TLE
%   (Note that TLEs provide mean-of-date data, not osculating)
M = deg2rad(325.0288); % Mean anomaly
n = 15.72125391*2*pi/86400; % Mean motion (rev/day) --> (rad/s)
a = gravity.GM^(1/3)*n^(-2/3); % Semimajor axis
e = 0.0006703; % Orbital eccentricity
I = deg2rad(51.6416); % Orbital inclination
w = deg2rad(130.5360); % Argument of periapsis
W = deg2rad(247.4627); % Right ascension of ascending node (RAAN)
tau = -M/n; % Last periapsis time with t0 = 0

%% ODE Propagation
% Solve for the true anomaly via Kepler's equation to obtain the true
% anomaly
[~, v] = anomalies(M, e);
% Convert TLE to modified equinoctial elements
p = a*(1 - e^2);
f = e*cos(w + W);
g = e*sin(w + W);
h = tan(I/2)*cos(W);
k = tan(I/2)*sin(W);
L = W + w + v;

tspan = [0, 50*(2*pi/n)]; % Time span (~x complete orbits)
x0 = [p, f, g, h, k, L]'; % Initial state of modified equinoctial elements
options = odeset('InitialStep', 15, 'RelTol', 1e-13, 'MaxStep', 30, 'Events', @odevents);

% Solve
sol = ode113(@odeval, tspan, x0, options, gravity);

% Turn the solution to be a time-history, where time flows downwards
sol.x = sol.x';
sol.y = sol.y';
sol.ie = sol.ie';
sol.xe = sol.xe';
sol.ye = sol.ye';

% Solution key
t = sol.x;
p = sol.y(:,1);
f = sol.y(:,2);
g = sol.y(:,3);
h = sol.y(:,4);
k = sol.y(:,5);
L = sol.y(:,6);

%% Plot
% Obtain inertial position
Reci = modeq2inertial(p, f, g, h, k, L);
% Plot over the Earth
satglobe4e, hold on
plot3(Reci(:,1), Reci(:,2), Reci(:,3), 'Color', '#A2142F')
plot3(Reci(end,1), Reci(end,2), Reci(end,3), 'w.', 'MarkerSize', 5), hold off
% Plot radius
figure, plot(t/86400, sqrt(Reci(:,1).^2+Reci(:,2).^2+Reci(:,3).^2)-6378.1363), axis tight
xlabel("Time $t$ [days]", 'interpreter', 'latex')
ylabel("Altitude $h$ [km]", 'interpreter', 'latex')