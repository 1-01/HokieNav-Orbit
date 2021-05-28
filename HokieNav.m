% HokieNav Orbital Simulation
% Matt Werner (m.werner@vt.edu) - April 8, 2021
clear, clc

%% Inputs
% Gravitational parameter (km3/s2)
gravity.GM = 398600.4418;

% Gravity model degree and order
[gravity.degree, gravity.order, ~, gravity.Cnm, gravity.Snm, ~, ~] ...
    = loadGravitationalCoefficients(10, 10, "tide-free");

% Initial time
[IC.t.JD, IC.t.UTC, IC.t.local] = defineInitialTime('20-May-2021 04:06:02.000', 'd-MMMM-yyyy HH:mm:ss.SSS', '-00:00');

TLE.M = deg2rad(89.9422); % Mean anomaly
TLE.n = 15.49024471*2*pi/86400; % Mean motion (rev/day) --> (rad/s)
TLE.a = gravity.GM^(1/3)*TLE.n^(-2/3); % Semimajor axis
TLE.e = 0.0003181; % Orbital eccentricity
TLE.I = deg2rad(51.6437); % Orbital inclination
TLE.w = deg2rad(20.5176); % Argument of periapsis
TLE.W = deg2rad(121.9403); % Right ascension of ascending node (RAAN)

%% ODE Propagation
% Solve for the true anomaly via Kepler's equation
[~, TLE.v] = anomalies(TLE.M, TLE.e);
% Convert TLE to modified equinoctial elements
[IC.MEE.p, IC.MEE.f, IC.MEE.g, IC.MEE.h, IC.MEE.k, IC.MEE.L] ...
    = orbital2modeq(TLE.a, TLE.e, TLE.I, TLE.w, TLE.W, TLE.v);

tspan = [0, 15.5*(2*pi/TLE.n)]; % Time span (~x complete orbits)
tspan = [0, 50405];
x0 = struct2array(IC.MEE)'; % Initial state of modified equinoctial elements
options = odeset('InitialStep', 15, 'RelTol', 1e-13, 'MaxStep', 30, 'Events', @odevents);

% Solve
traj = struct('t', [], 'MEE', [], 'events', []);
[traj.t, MEE, traj.events.t, traj.events.MEE, traj.events.ID] ...
    = ode113(@odeval, tspan, x0, options, gravity, IC);
traj.MEE.p = MEE(:,1); traj.MEE.f = MEE(:,2); traj.MEE.g = MEE(:,3);
traj.MEE.h = MEE(:,4); traj.MEE.k = MEE(:,5); traj.MEE.L = MEE(:,6);
clear MEE

%% Dynamics Reconstruction
for k = numel(traj.t):-1:1
    [~, ...
     traj.dMEEdt.p(k,1), ... 1
     traj.dMEEdt.f(k,1), ... 2
     traj.dMEEdt.g(k,1), ... 3
     traj.dMEEdt.h(k,1), ... 4
     traj.dMEEdt.k(k,1), ... 5
     traj.dMEEdt.L(k,1), ... 6
     traj.accel.grav(k,:), ... 7
     traj.accel.drag(k,:), ... 8
     traj.atmos.density(k,1), ... 9
     traj.timeSystems.GMST(k,1), ... 10
     traj.position.geodetic.longitude(k,1), ... 11
     traj.position.geodetic.latitude(k,1), ... 12
     traj.position.geodetic.altitude(k,1), ... 13
     traj.position.geocentric.radius(k,1), ... 14
     traj.position.geocentric.colatitude(k,1), ... 15
     traj.position.geocentric.longitude(k,1), ... 16
     traj.position.ECF.x(k,1), ... 17
     traj.position.ECF.y(k,1), ... 18
     traj.position.ECF.z(k,1), ... 19
     traj.position.ECI.X(k,1), ... 20
     traj.position.ECI.Y(k,1), ... 21
     traj.position.ECI.Z(k,1)] ... 22
        = odeval(traj.t(k), ...
              [traj.MEE.p(k), traj.MEE.f(k), traj.MEE.g(k), ...
               traj.MEE.h(k), traj.MEE.k(k), traj.MEE.L(k)], ...
              gravity, IC);
end

%% Additional Derived Information
% Recover the classical orbital elements
[traj.COE.a, traj.COE.e, traj.COE.I, traj.COE.w, traj.COE.W, traj.COE.v] ...
     = modeq2orbital(traj.MEE.p, traj.MEE.f, traj.MEE.g, traj.MEE.h, traj.MEE.k, traj.MEE.L);

% Perform coordinate transformations to express known quantities in
% different bases. Of particular interest is obtaining the Earth-fixed
% and Earth-inertial representations of the perturbing accelerations
Tsphr_orbsphr = [1,0,0;0,0,-1;0,1,0];
for k = numel(traj.t):-1:1
    Tcart_sphr = sphr2ijk(traj.position.geocentric.longitude(k), ...
                          traj.position.geocentric.colatitude(k));
    traj.accel.gravcart(k,:) = (Tcart_sphr*Tsphr_orbsphr*traj.accel.grav(k,:)')';
    traj.accel.dragcart(k,:) = (Tcart_sphr*Tsphr_orbsphr*traj.accel.grav(k,:)')';
    traj.accel.acart(k,:) = traj.accel.gravcart(k,:) + traj.accel.dragcart(k,:);
end
% Clear all temporary variables
clear k Tcart_sphr Tsphr_orbsphr tspan x0

%% Plot
% Plot over the Earth
satglobe4e, hold on
plot3(traj.position.ECF.x, ...
      traj.position.ECF.y, ...
      traj.position.ECF.z, ...
      'Color', '#A2142F')
plot3(traj.position.ECF.x(end), ...
      traj.position.ECF.y(end), ...
      traj.position.ECF.z(end), ...
      'w.', 'MarkerSize', 5), hold off
% Plot GPS height
figure, plot(traj.t/86400, traj.position.geodetic.altitude), axis tight
xlabel("Time $t$ [days]", 'interpreter', 'latex')
ylabel("Altitude $h$ [km]", 'interpreter', 'latex')