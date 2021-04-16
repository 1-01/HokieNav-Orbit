% HokieNav Orbital Simulation
% Matt Werner (m.werner@vt.edu) - April 8, 2021
clear, clc

%% Inputs
% Gravitational parameter
GM = getGM(1000);

% Initial conditions read from TLE
%   (Note that TLEs provide mean-of-date data, not osculating)
M = deg2rad(325.0288); % Mean anomaly
n = 15.72125391*2*pi/86400; % Mean motion (rev/day) --> (rad/s)
a = GM^(1/3)*n^(-2/3); % Semimajor axis
e = 0.0006703; % Orbital eccentricity
I = deg2rad(51.6416); % Orbital inclination
w = deg2rad(130.5360); % Argument of periapsis
W = deg2rad(247.4627); % Right ascension of ascending node (RAAN)
tau = -M/n; % Last periapsis time with t0 = 0

%% ODE Propagation
tspan = [0, 1*(2*pi/n)]; % Time span (1 complete orbit)
x0 = [a, e, I, w, W, tau]'; % Initial state of classical orbital elements
options = odeset('RelTol', 1e-7, 'MaxStep', 30, 'Events', @odevents); % Options

% Solve
sol = ode113(@odeval, tspan, x0, options);

% Turn the solution to be a time-history, where time flows downwards
sol.x = sol.x';
sol.y = sol.y';
sol.ie = sol.ie';
sol.xe = sol.xe';
sol.ye = sol.ye';

% Solution key
t = sol.x;
a = sol.y(:,1);
e = sol.y(:,2);
I = sol.y(:,3);
w = sol.y(:,4);
W = sol.y(:,5);
tau = sol.y(:,6);

%% Plot
M = meanMotion(GM, a).*(t - tau);
[~, f] = anomalies(M, e);
% Compute inertial position
Reci = orbital2inertial([a, e, I, w, W, f]);

% Plot over the Earth
satglobe4e, hold on, plot3(Reci(:,1), Reci(:,2), Reci(:,3), 'Color', '#A2142F'), plot3(Reci(end,1), Reci(end,2), Reci(end,3), '.', 'MarkerSize', 5, 'Color', 'w')