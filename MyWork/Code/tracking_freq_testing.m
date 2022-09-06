%% Housekeeping
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

%% Input data

% te in Days
% R, SMA in km
% V in km/s
% Acceleration in km/s^2
% E and Chi in deg
% X Y Z in km
% Temp in K

[te, R, SMA, V, Ax, Ay, Az, E1, E2, E3, E1_dot, E2_dot, E3_dot, X, Y, Z, chi2, chi3, Temp, Tmax, th, Ve, Amax, Rmin] = ImportGMATData(1);
te_s = te*24*60*60; % Time elapsed in seconds

%% Closest approach timing
closest_approach_index = find(R==Rmin);
R_closest = R(closest_approach_index);
te_closest = te(closest_approach_index);

%% Max tracking - need to be able to track until
R_maxtrack = 8*AU;
maxtrack_index = find(R<R_maxtrack, 1, 'last');
te_maxtrack = floor(te(maxtrack_index));

%% Consumption parameters
passive_consumption = 0.4; %Watts
% During transmission (over estimating length of transmission to 1h)
max_consumption = 5; % Watts

% Assuming datarate of 512 bps with transmission time of 30 minutes or 1h
% minute giving a file size of 0.9 mb or 1.8 mb
link_time_min = 30;
link_time_hr = link_time_min/60;

% Assume voltage of amplifier is 3.3V
link_energy_usage = max_consumption*link_time_hr;

% energy consumption per day
link_energy_day = max_consumption*link_time_hr + passive_consumption*(day-link_time_hr);
passive_energy_day = passive_consumption*day;

for j = 1:3
    %% Energy consumption based on tracking rates
    tracking_rate = j; % Pick tracking rate here

    weekly = 1:7:te_maxtrack;
    monthly = 1:28:te_maxtrack;

    % initialise variables
    t_tracking = 1:1:te_maxtrack;
    energy_con = ones(te_maxtrack,1);
    switch tracking_rate
    % Daily - 1
        case 1
            energy_con = energy_con*link_energy_day;
    % Weekly - 2
        case 2
            energy_con = energy_con*passive_energy_day;
            energy_con(weekly) = link_energy_day;
    % Monthly (every 28 days) - 3
        case 3
            energy_con = energy_con*passive_energy_day;
            energy_con(monthly) = link_energy_day;
    end

    %% Cumulative energy consumption
    energy_con_cum = zeros(te_maxtrack,1);
    energy_con_cum(1) = link_energy_day;
    for i = 2:te_maxtrack
        energy_con_cum(i) = energy_con_cum(i-1) + energy_con(i);
    end
    tot_energy_con = energy_con_cum(end);

    %% Plotting consumption
    figure(1)
    hold on
    plot(t_tracking/yr, energy_con,'LineWidth',0.01)
    grid on
    box on
    title('Energy Consumption each day')
    xlabel('Years elapsed [yrs]')
    ylabel('Energy conusmed each day [W*hrs]')
    legend('Daily','Weekly','Monthly','Location','Northwest')
    axis([0 8000/yr 8 13])
    hold off

    figure(2)
    hold on
    plot(t_tracking/yr, energy_con_cum)
    grid on
    box on
    title('Cumulative Energy Consumption')
    xlabel('Years elapsed [yrs]')
    ylabel('Total energy conusmed [W*hrs]')
    legend('Daily','Weekly','Monthly','Location','Southeast')
    axis([0 8000/yr 0 9*10^4])
    hold off
    
    clear energy_con energy_con_cum 
end