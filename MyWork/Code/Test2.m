%% houskeeping
clear all
close all
clc

%% Constants
% AU to km
AU = 149597870.7; % [km]
% Speed of Light
c = 299792458; % [m/s]
% Years to days
yr = 365.2422; % [days]
% Hours in a day
day = 24;

%% Inputs
% Dimension of Sail
side = sqrt(28.6);

% Rotation of body about spin axis
angle_spinaxis = 180; % Orientation about spin axis (clockwise)
theta = 360-angle_spinaxis; % Anti-clockwise value

% Angle between Sail and Sun (depends on TA + attitude propagation)
alpha = 10; 

%% SRP over length (1D Case)
%

dist = abs(side*cosd(theta) - side*sind(theta));

local_coords = linspace(-dist/2,dist/2)*cosd(80);

global_coords = AU*0.13*1000 + local_coords;

SRP_local = SolarIntensity(global_coords)/c; % SRP over length of sail

range = linspace(0,8*AU);

SRP_all = SolarIntensity(range*1000)/c; 

%% Plot
figure(1)
plot(local_coords, SRP_local*10^(6))
xlabel('Local X [m]')
ylabel('SRP [\muPa]')


%{
figure(2)
plot(range/AU, SRP_all*10^(6))
xlabel('Range [AU]')
ylabel('SRP [\muPa]')
%}