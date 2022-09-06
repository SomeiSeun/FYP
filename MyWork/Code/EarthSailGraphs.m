% Code to generate radius from Earth Graph
%% housekeeping

clear all
close all
clc

%% Constants
% AU to km
AU = 149597870.7; % [km]
% Years to days
yr = 365.2422; % [days]
% hours in a day
day = 24;
% Solar constants
I0 = 1365.4; % Solar intensity [W/m^2]
Cr = [0.91,0.16]; % Reflectivity coefficient vector for the two
% Solar gravitational constant [m^3/s^2]
mu = 1.32712440018e20;
% sail surfaces: [reflective,emissive] [-]
SBc = 5.670374419e-8; % Stefan-Boltzmann constant [W/(m^2 K^4)]
% Perihelion radius set for simulations [AU]
mass = [2.27,3,5,10,1.5]; % initial yaw angle

%% Import data
[te, R, SMA, V, Ax, Ay, Az, E1, E2, E3, E1_dot, E2_dot, E3_dot, X, Y, Z, chi2, chi3, G, Temp, Tmax, th, Ve, Amax, Rmin, R_earth] = ImportGMATData(1);

%% Calculating temperature

% Max tracking - need to be able to track until
R_maxtrack = 8*AU;
maxtrack_index = find(R<R_maxtrack, 1, 'last');
te_maxtrack = floor(te(maxtrack_index));
t_tracking = 1:1:te_maxtrack; t_tracking = t_tracking';

rho_r = 0.91; % Reflectivitiy of reflective surface
rho_e = 0.16; % Reflectivity of emissive surface
SBc = 5.670374419e-8; % Steffan boltzmann constant
R_tracking = interp1q(te,R,t_tracking);
chi2_daily = interp1q(te,chi2,t_tracking);
chi3_daily = interp1q(te,chi3,t_tracking);
Temp_sail_daily = Temp_sail(R_tracking*1000,chi2_daily,chi3_daily,rho_r,rho_e);


%% Plotting

limit_21yr = length(find((te < 21*yr)));

fig1 = figure(1);
plot(te/yr, R_earth/AU)
title('Sailcraft-Earth distance over mission')
xlabel('Time elapsed [yrs]')
xlim([0 te(limit_21yr)/yr])
ylabel('Distance from Earth [AU]')
box on
grid on
grid minor

fig1.Units = 'inches';
fig1.Position(3) = 4.5;
fig1.Position(4) = 3.2;
set(fig1.Children, 'FontName', 'Arial', 'FontSize', 11);
%print('EarthSailRadius', '-depsc')

fig2 = figure(2);
plot(te/yr, R/AU)
title('Sailcraft-Sun distance over mission')
xlabel('Time elapsed [yrs]')
ylabel('Distance from Sun [AU]')
box on
grid on
grid minor
ylim([0 123])
xlim([0 55])

fig2.Units = 'inches';
fig2.Position(3) = 4.5;
fig2.Position(4) = 3.2;
set(fig2.Children, 'FontName', 'Arial', 'FontSize', 11);
print('SunSailRadius', '-depsc')

fig3 = figure(3);
plot(te/yr, V)
title('Sailcraft velocity over mission')
xlabel('Time elapsed [yrs]')
ylabel('Velocity of Sailcraft [kms^{-1}]')
box on
grid on
grid minor
xlim([0 55])

fig3.Units = 'inches';
fig3.Position(3) = 4.5;
fig3.Position(4) = 3.2;
set(fig3.Children, 'FontName', 'Arial', 'FontSize', 11);
print('SunSailVelocity', '-depsc')

fig4 = figure(4);
hold on
plot(t_tracking/yr,Temp_sail_daily)
grid on
grid minor
box on
title('Temperature of sailcraft over mission')
xlabel('Time elapsed [yrs]')
ylabel('Sailcraft temperature [K]')
xlim([0 21])
hold off

fig4.Units = 'inches';
fig4.Position(3) = 4.5;
fig4.Position(4) = 3.2;
set(fig4.Children, 'FontName', 'Arial', 'FontSize', 11);
print('TempTime', '-depsc')
