% HokieNav Orbital Simulation
% Matt Werner (m.werner@vt.edu) - April 8, 2021
clear, clc, HokieNavPath

%% Inputs
diary out/log
disp("HokieNav Orbital Simulation")
disp(strcat("Start time: ", datestr(datetime('now'))))

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
[~, TLE.f] = anomalies(TLE.M, TLE.e);
% Convert TLE to modified equinoctial elements
[IC.MEE.p, IC.MEE.f, IC.MEE.g, IC.MEE.h, IC.MEE.k, IC.MEE.L] ...
    = orbital2modeq(TLE.a, TLE.e, TLE.I, TLE.w, TLE.W, TLE.f);
% Convert the initial conditions to the (ECI) Cartesian state
[IC.x.r, IC.x.v] = modeq2inertial(IC.MEE.p, IC.MEE.f, IC.MEE.g, IC.MEE.h, IC.MEE.k, IC.MEE.L, gravity.GM);

% Time span
tspan = [0, 15.5*(2*pi/TLE.n)]; % (~x complete orbits)
tspan = [0, 50405];
tspan = [0, 600];

% Initial condition
x0 = struct2array(IC.x)'; % Initial state of modified equinoctial elements

% Specify options for the ODE solver
options = odeset('InitialStep', 15, ...
                 'RelTol', 1e-13, ...
                 'AbsTol', 1e-14, ...
                 'MaxStep', 0.01, ...
                 'Events', @odevents);

% Solve
traj = struct('t', [], 'x', [], 'dxdt', [], 'events', []);
[traj.t, x, traj.events.t, traj.events.x, traj.events.ID] ...
    = ode113(@odeval, tspan, x0, options, gravity, IC);

% Distribute the state
traj.x.PX = x(:,1); traj.x.PY = x(:,2); traj.x.PZ = x(:,3);
traj.x.VX = x(:,4); traj.x.VY = x(:,5); traj.x.VZ = x(:,6);
clear x

%% Dynamics Reconstruction
for k = numel(traj.t):-1:1
    [~, ...
     traj.dxdt.VX(k,1), ... 1
     traj.dxdt.VY(k,1), ... 2
     traj.dxdt.VZ(k,1), ... 3
     traj.dxdt.AX(k,1), ... 4
     traj.dxdt.AY(k,1), ... 5
     traj.dxdt.AZ(k,1), ... 6
     traj.accelerations.ECI.gravity(k,:), ... 7
     traj.accelerations.ECI.drag(k,:), ... 8
     traj.atmosphere.density(k,1), ... 9
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
                  [traj.x.PX(k); traj.x.PY(k); traj.x.PZ(k);  ...
                   traj.x.VX(k); traj.x.VY(k); traj.x.VZ(k)], ...
                  gravity, IC);
end

%% Additional Derived Information

% Perform coordinate transformations to express known quantities in
% different bases. Of particular interest is obtaining the Earth-fixed
% and spherical representations of the perturbing accelerations
for k = numel(traj.t):-1:1
    % Recover the classical orbital elements
    [traj.COE.a(k,1), traj.COE.e(k,1), traj.COE.I(k,1), ...
     traj.COE.w(k,1), traj.COE.W(k,1), traj.COE.f(k,1)] ...
      = inertial2orbital([traj.x.PX(k); traj.x.PY(k); traj.x.PZ(k)], ...
                         [traj.x.VX(k); traj.x.VY(k); traj.x.VZ(k)], gravity.GM);
    
    % Coordinate transformations
    TECF_ECI = R3(traj.timeSystems.GMST(k));
    TSPH_ECF = sphr2ijk(traj.position.geocentric.longitude(k), ...
                          traj.position.geocentric.colatitude(k))';
    traj.accelerations.ECF.gravity(k,:) = (TECF_ECI*traj.accelerations.ECI.gravity(k,:)')';
    traj.accelerations.SPH.gravity(k,:) = (TSPH_ECF*traj.accelerations.ECF.gravity(k,:)')';
    traj.accelerations.ECF.drag(k,:) = (TECF_ECI*traj.accelerations.ECI.drag(k,:)')';
    traj.accelerations.SPH.drag(k,:) = (TSPH_ECF*traj.accelerations.ECF.drag(k,:)')';
end
% Clear all temporary variables
clear k TSPH_ECF TECF_ECI tspan x0
diary off

%% Plot
% Plot over the Earth
plotTrajectory(traj.position.ECI.X, ...
               traj.position.ECI.Y, ...
               traj.position.ECI.Z)

% Plot ground track
figure
groundTrack(traj.position.geodetic.longitude, ...
            traj.position.geodetic.latitude)