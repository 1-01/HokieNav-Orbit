GMAT = readtable("C:\Users\Matt\AppData\Local\GMAT\R2018a\output\GMAT_SIM_ISS_2021_May_20_04_06_02_UTC_(Day_140).txt");

% Interpolate Matlab results to plot residual
MATLAB.t = interp1(traj.t/86400, traj.t/86400, GMAT.HokieNav_ElapsedDays);
MATLAB.PX = interp1(traj.t/86400, traj.x.PX, MATLAB.t);
MATLAB.PY = interp1(traj.t/86400, traj.x.PY, MATLAB.t);
MATLAB.PZ = interp1(traj.t/86400, traj.x.PZ, MATLAB.t);
MATLAB.VX = interp1(traj.t/86400, traj.x.VX, MATLAB.t);
MATLAB.VY = interp1(traj.t/86400, traj.x.VY, MATLAB.t);
MATLAB.VZ = interp1(traj.t/86400, traj.x.VZ, MATLAB.t);

% Plot residuals of radius and velocity
subplot(2,1,1)
plot(MATLAB.t, sqrt((MATLAB.PX - GMAT.HokieNav_EarthMJ2000Eq_X).^2 + ...
                    (MATLAB.PY - GMAT.HokieNav_EarthMJ2000Eq_Y).^2 + ...
                    (MATLAB.PZ - GMAT.HokieNav_EarthMJ2000Eq_Z).^2))
hold on
% % plot(MATLAB.t, sqrt(cumtrapz(MATLAB.t, abs(MATLAB.PX - GMAT.HokieNav_EarthMJ2000Eq_X).^2)))
% % plot(MATLAB.t, sqrt(cumtrapz(MATLAB.t, abs(MATLAB.PY - GMAT.HokieNav_EarthMJ2000Eq_Y).^2)))
% % plot(MATLAB.t, sqrt(cumtrapz(MATLAB.t, abs(MATLAB.PZ - GMAT.HokieNav_EarthMJ2000Eq_Z).^2)))
% % plot(MATLAB.t, sqrt(cumtrapz(MATLAB.t, sqrt((MATLAB.PX - GMAT.HokieNav_EarthMJ2000Eq_X).^2 + (MATLAB.PY - GMAT.HokieNav_EarthMJ2000Eq_Y).^2 + (MATLAB.PZ - GMAT.HokieNav_EarthMJ2000Eq_Z).^2).^2)))
% % % plot(MATLAB.t, abs(sqrt(MATLAB.PX.^2 + MATLAB.PY.^2 + MATLAB.PZ.^2) - sqrt(GMAT.HokieNav_EarthMJ2000Eq_X.^2 + GMAT.HokieNav_EarthMJ2000Eq_Y.^2 + GMAT.HokieNav_EarthMJ2000Eq_Z.^2)))
ylabel("Position Residual [km]", 'interpreter', 'latex')
hold off

subplot(2,1,2)
plot(MATLAB.t, sqrt((MATLAB.VX - GMAT.HokieNav_EarthMJ2000Eq_VX).^2 + ...
                    (MATLAB.VY - GMAT.HokieNav_EarthMJ2000Eq_VY).^2 + ...
                    (MATLAB.VZ - GMAT.HokieNav_EarthMJ2000Eq_VZ).^2))
% % hold on
% % plot(MATLAB.t, sqrt(cumtrapz(MATLAB.t, abs(MATLAB.VX - GMAT.HokieNav_EarthMJ2000Eq_VX).^2)))
% % plot(MATLAB.t, sqrt(cumtrapz(MATLAB.t, abs(MATLAB.VY - GMAT.HokieNav_EarthMJ2000Eq_VY).^2)))
% % plot(MATLAB.t, sqrt(cumtrapz(MATLAB.t, abs(MATLAB.VZ - GMAT.HokieNav_EarthMJ2000Eq_VZ).^2)))
% % plot(MATLAB.t, sqrt(cumtrapz(MATLAB.t, sqrt((MATLAB.VX - GMAT.HokieNav_EarthMJ2000Eq_VX).^2 + (MATLAB.VY - GMAT.HokieNav_EarthMJ2000Eq_VY).^2 + (MATLAB.VZ - GMAT.HokieNav_EarthMJ2000Eq_VZ).^2).^2)))
xlabel("Time $t$ [days]", 'interpreter', 'latex')
ylabel("Velocity Residual [km/s]", 'interpreter', 'latex')
hold off

% Plot MJ2000 results from GMAT
figure
subplot(2,1,1)
plot(GMAT.HokieNav_ElapsedDays, [GMAT.HokieNav_EarthMJ2000Eq_X,  GMAT.HokieNav_EarthMJ2000Eq_Y,  GMAT.HokieNav_EarthMJ2000Eq_Z])
ylabel("MJ2000 Position", 'interpreter', 'latex')
hold on
subplot(2,1,2)
plot(GMAT.HokieNav_ElapsedDays, [GMAT.HokieNav_EarthMJ2000Eq_VX, GMAT.HokieNav_EarthMJ2000Eq_VY, GMAT.HokieNav_EarthMJ2000Eq_VZ])
xlabel("Time $t$ [days]", 'interpreter', 'latex')
ylabel("MJ2000 Velocity", 'interpreter', 'latex')
hold on

% Plot J2000 results from Matlab
subplot(2,1,1)
plot(traj.t/86400, [traj.x.PX, traj.x.PY, traj.x.PZ])
ylabel("J2000 Position", 'interpreter', 'latex')
hold off
subplot(2,1,2)
plot(traj.t/86400, [traj.x.VX, traj.x.VY, traj.x.VZ])
xlabel("Time $t$ [days]", 'interpreter', 'latex')
ylabel("MJ2000 Velocity", 'interpreter', 'latex')
hold off
% % Position...
% % xlim: 0.159307565574916         0.159309985297721
% % ylim: 3033.94995119241          3044.42224454063