function traj = propagateTrajectory(TLEfile)
% 
% Matt Werner (m.werner@vt.edu) - April 8, 2021
% 
% Simulate the trajectory of a CubeSat using the last known two-line
% element set (TLE) before the satellite undergoes GNSS communication
% blackout. The nominal length of time that the satellite will be without
% communications is ~10 minutes (~600 seconds). 
% 
% During this period of time, the satellite may transit over a measurable
% portion of Earth by undergoing standard Keplerian motion as well as
% encountering perturbations (gravity, drag, SRP, etc.).
% 
%    Inputs:
% 
%           TLEfile - File name of the two-line element set containing the
%                     most recent information about the satellite's
%                     trajectory before GNSS blackout. This file contains
%                     standard meta data (satellite name, etc.) as well as
%                     initial conditions (time, position, and velocity) in
%                     the form of the Keplerian orbital elements. Units of
%                     quantities within the TLE are standard. That is, 
% 
%    Outputs:
% 
%              traj - Time-history of the trajectory since the epoch
%                     provided by the given TLE file.
% 

%% Gravity parameters
% Gravitational parameter (km3/s2)
gravity.GM = 398600.4418;

% Gravity model degree and order
[gravity.degree, gravity.order, ~, gravity.Cnm, gravity.Snm, ~, ~] ...
    = loadGravitationalCoefficients(10, 10, 'tide-free');

%% Initial conditions
% Read the initial time and the initial state (orbital elements) from TLE
[IC.t, IC.COE] = tle(TLEfile, 'rad');

% Convert the orbital elements to modified equinoctial elements and then to
% the (ECI) Cartesian state
IC.MEE = orbital2modeq(IC.COE);
IC.ECI = modeq2inertial(IC.MEE, gravity.GM);
%% ODE Propagation
% Specify options for the ODE solver
options = odeset('InitialStep', 1, ...
                 'RelTol', 1e-13, ...
                 'AbsTol', 1e-14, ...
                 'MaxStep', inf, ...
                 'Events', @odevents);
% Time span (Note that the actual timespan we want is 10 minutes from
% *now*, but the latest TLE update might be some time old so that the
% effective integration time is more than 10 minutes.)
tspan = [0, 600];
% Initial condition for the state x
x0 = struct2array(IC.ECI)';

% Solve
traj = struct('IC', [], 't', [], 'x', [], 'dxdt', [], 'events', []);
[traj.t, x, traj.events.t, traj.events.x, traj.events.ID] ...
    = ode113(@odeval, tspan, x0, options, gravity, IC);

% Distribute the state
traj.IC = IC;
traj.x.PX = x(:,1); traj.x.PY = x(:,2); traj.x.PZ = x(:,3);
traj.x.VX = x(:,4); traj.x.VY = x(:,5); traj.x.VZ = x(:,6);
clear x

%% Dynamics Reconstruction
% Recover information that was previously calculated in obtaining the
% solution
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
    tmpECI.X = traj.x.PX(k);
    tmpECI.Y = traj.x.PY(k);
    tmpECI.Z = traj.x.PZ(k);
    tmpECI.VX = traj.x.VX(k);
    tmpECI.VY = traj.x.VY(k);
    tmpECI.VZ = traj.x.VZ(k);
    tmpCOE = inertial2orbital(tmpECI, gravity.GM);
    traj.COE.a(k,1) = tmpCOE.a;
    traj.COE.e(k,1) = tmpCOE.e;
    traj.COE.I(k,1) = tmpCOE.I;
    traj.COE.W(k,1) = tmpCOE.W;
    traj.COE.w(k,1) = tmpCOE.w;
    traj.COE.f(k,1) = tmpCOE.f;
    
    % Coordinate transformations
    TECF_ECI = R3(traj.timeSystems.GMST(k));
    TSPH_ECF = sphr2ijk(traj.position.geocentric.longitude(k), traj.position.geocentric.colatitude(k))';
    traj.accelerations.ECF.gravity(k,:) = (TECF_ECI*traj.accelerations.ECI.gravity(k,:)')';
    traj.accelerations.SPH.gravity(k,:) = (TSPH_ECF*traj.accelerations.ECF.gravity(k,:)')';
    traj.accelerations.ECF.drag(k,:) = (TECF_ECI*traj.accelerations.ECI.drag(k,:)')';
    traj.accelerations.SPH.drag(k,:) = (TSPH_ECF*traj.accelerations.ECF.drag(k,:)')';
end

%% Plot
% Plot over the Earth
plotTrajectory(traj.position.ECI.X, ...
               traj.position.ECI.Y, ...
               traj.position.ECI.Z)

% Plot ground track
figure
groundTrack(traj.position.geodetic.longitude, ...
            traj.position.geodetic.latitude)