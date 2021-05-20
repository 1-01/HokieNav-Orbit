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
TLE.M = deg2rad(325.0288); % Mean anomaly
TLE.n = 15.72125391*2*pi/86400; % Mean motion (rev/day) --> (rad/s)
TLE.a = gravity.GM^(1/3)*TLE.n^(-2/3); % Semimajor axis
TLE.e = 0.0006703; % Orbital eccentricity
TLE.I = deg2rad(51.6416); % Orbital inclination
TLE.w = deg2rad(130.5360); % Argument of periapsis
TLE.W = deg2rad(247.4627); % Right ascension of ascending node (RAAN)
% tau = -M/n; % Last periapsis time with t0 = 0

TLE.M = deg2rad(89.9422); % Mean anomaly
TLE.n = 15.49024471284232*2*pi/86400; % Mean motion (rev/day) --> (rad/s)
TLE.a = gravity.GM^(1/3)*TLE.n^(-2/3); % Semimajor axis
TLE.e = 0.0003181; % Orbital eccentricity
TLE.I = deg2rad(51.6437); % Orbital inclination
TLE.w = deg2rad(20.5176); % Argument of periapsis
TLE.W = deg2rad(121.9403); % Right ascension of ascending node (RAAN)
% tau = -M/n; % Last periapsis time with t0 = 0

%% ODE Propagation
% Solve for the true anomaly via Kepler's equation to obtain the true
% anomaly
[~, TLE.v] = anomalies(TLE.M, TLE.e);
% Convert TLE to modified equinoctial elements
IC.p = TLE.a*(1 - TLE.e^2);
IC.f = TLE.e*cos(TLE.w + TLE.W);
IC.g = TLE.e*sin(TLE.w + TLE.W);
IC.h = tan(TLE.I/2)*cos(TLE.W);
IC.k = tan(TLE.I/2)*sin(TLE.W);
IC.L = TLE.W + TLE.w + TLE.v;
 
tspan = [0, 15.5*(2*pi/TLE.n)]; % Time span (~x complete orbits)
x0 = struct2array(IC)'; % Initial state of modified equinoctial elements
options = odeset('InitialStep', 15, 'RelTol', 1e-13, 'MaxStep', 30, 'Events', @odevents);

% Solve
[sol.t, sol.x, sol.te, sol.xe, sol.ie] = ode113(@odeval, tspan, x0, options, gravity);

%% Plot
% Obtain inertial position
[a, e, I, w, W, v] ...
     = modeq2orbital(sol.x(:,1), sol.x(:,2), sol.x(:,3), sol.x(:,4), sol.x(:,5), sol.x(:,6));
Reci = modeq2ecf(sol.x(:,1), sol.x(:,2), sol.x(:,3), sol.x(:,4), sol.x(:,5), sol.x(:,6));
% Plot over the Earth
satglobe4e, hold on
plot3(Reci(:,1), Reci(:,2), Reci(:,3), 'Color', '#A2142F')
plot3(Reci(end,1), Reci(end,2), Reci(end,3), 'w.', 'MarkerSize', 5), hold off
% Plot radius
figure, plot(sol.t/86400, sqrt(Reci(:,1).^2 + Reci(:,2).^2 + Reci(:,3).^2) - 6378.1363), axis tight
xlabel("Time $t$ [days]", 'interpreter', 'latex')
ylabel("Altitude $h$ [km]", 'interpreter', 'latex')