% Code to compute power generation vs power consumption rates over course
% of entire journey
% Code by Tanmay | 220424

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


%% Input data

% te in Days
% R, SMA in km
% V in km/s
% Acceleration in km/s^2
% E and Chi in deg
% X Y Z in km
% Temp in K

[te, R, SMA, V, Ax, Ay, Az, E1, E2, E3, E1_dot, E2_dot, E3_dot, X, Y, Z, chi2, chi3, G, Temp, Tmax, th, Ve, Amax, Rmin, R_earth] = ImportGMATData(1);
TA = G(:,2); % True anomaly in degrees

te_s = te*24*60*60; % Time elapsed in seconds
n = length(te);

% chi3 to be used for true anomaly

%% Power Consumption
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% POWER CONSUMPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Closest approach timing
closest_approach_index = find(R==Rmin);
R_closest = R(closest_approach_index);
te_closest = te(closest_approach_index);

%% Max tracking - need to be able to track until
R_maxtrack = 8*AU;
maxtrack_index = find(R<R_maxtrack, 1, 'last');
te_maxtrack = floor(te(maxtrack_index));

%% Consumption parameters
passive_consumption = 0.2+0.25; %Watts
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

% initialise variables
t_tracking = 1:1:te_maxtrack; t_tracking = t_tracking';
weekly = 1:7:te_maxtrack;
monthly = 1:28:te_maxtrack;
energy_con = ones(te_maxtrack,3);
for j = 1:3
    %% Energy consumption based on tracking rates
    tracking_rate = j; % Pick tracking rate here
    
    switch tracking_rate
    % Daily - 1
        case 1
            energy_con(:,j) = energy_con(:,j)*link_energy_day;
    % Weekly - 2
        case 2
            energy_con(:,j) = energy_con(:,j)*passive_energy_day;
            energy_con(weekly,j) = link_energy_day;
    % Monthly (every 28 days) - 3
        case 3
            energy_con(:,j) = energy_con(:,j)*passive_energy_day;
            energy_con(monthly,j) = link_energy_day;
    end
end
%% Cumulative energy consumption
energy_con_cum = cumtrapz(t_tracking,energy_con);
tot_energy_con = trapz(t_tracking,energy_con);
%% Plotting consumption
fig1 = figure(1);
hold on
plot(t_tracking/yr, energy_con_cum(:,1),'-')
plot(t_tracking/yr, energy_con_cum(:,2),'--')
plot(t_tracking/yr, energy_con_cum(:,3),'-.k')
grid on
grid minor
box on
title('Cumulative Energy Consumption')
xlabel('Time elapsed [yrs]')
ylabel('Total energy conusmed [Wh]')
axis([0 21 0 10^5])
hold off
legend('Daily','Weekly','Monthly','Location','Southeast')

fig1.Units = 'inches';
fig1.Position(3) = 6;
fig1.Position(4) = 3;
set(fig1.Children, 'FontName', 'Arial', 'FontSize', 11);
print('PowerConsumption', '-depsc')

%% Only battery power
LiCFx_energydensity = 300; % Wh/kg
batteryonly_mass = tot_energy_con/LiCFx_energydensity; % in kg

% These values are absurdly high and would require a SmallSat or larger,
% and would completely invalidate the use of Solar Sails as the mass would
% be too high to impart sufficient impulse to escape the solar system
% within the required time constraints.

%% Baseline battery capacity

% A battery capacity needs to be decided upon 
batteryonly_mass = 1; % in kg
battery_capacity = batteryonly_mass*LiCFx_energydensity;


%% Power Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% POWER GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CubeSat Parameters
% Area (3U CubeSat)
dim_width = 0.1; % in m
dim_height = 0.3; % in m
area_sp = dim_width*dim_height*sqrt(2); % of each side panel in m^2

%% Solar Panels on side of CubeSat only

% Solar intensity at every calculated point
I = SolarIntensity(R*10^3);

% Will use Multi-Junction Solar Arrays for high efficiency making effective
% use of limited area avaliable
eff_mj = 0.22;

% Degradation rate of 3.75%
degradation_sp = 3.75/100;
L_d_sp = (1 - degradation_sp).^(te/365.2422); % Lifetime degredation

%Power generated by solar array
gen_sp = eff_mj.*I.*area_sp.*L_d_sp;

% Can generate power between 261.95 -> 98.05 degrees, else turned off 
solar_on1_sp = TA > 261.95;
solar_on2_sp = TA < 98.05;
solar_on_sp = solar_on1_sp + solar_on2_sp;

% When operational
gen_operational_sp = gen_sp.*solar_on_sp.*abs(cosd(TA));
energy_gen_cum_sp = cumtrapz(te,gen_operational_sp);

% Use figure 2 to showcase that power generated by solar panels on the
% doors is 1-2 orders of magnitude lower than needed for the total power
% consumption of the mission. As such, Thin Film Solar Cells must be
% utilised.

fig2 = figure(2);
hold on
plot(te/yr, energy_gen_cum_sp)
grid on
grid minor
box on
title('Cumulative Energy Generation - Solar Panels')
xlabel('Time elapsed [yrs]')
ylabel('Total energy generated [Wh]')
xlim([0 21])
hold off

fig2.Units = 'inches';
fig2.Position(3) = 6;
fig2.Position(4) = 3;
set(fig2.Children, 'FontName', 'Arial', 'FontSize', 11);
%print('SolarPanelGen', '-depsc')

%% Thin film solar cells in addition to solar panels

% CIGS on PI Substrate considered
CIGS_PI_mass_area = 41/1000; %kg/m^2

quad1 = TA < 90;
quad4 = TA > 270;
quad14 = quad1 + quad4; % Faces sun from 270-90 deg
quad23 = quad14 == 0; % Faces sun from 90-270 deg


area_sc1 = 0.1*quad14; % Faces sun from 270-90 deg
area_sc2 = 14.9704*quad23; % Faces sun from 90-270 deg 
area_sc = area_sc1 + area_sc2; % m^2
mass_sc = (0.1 + 14.9704)*CIGS_PI_mass_area;
% Using efficiency of 14.1% considering Na post-treatement method
eff_CIGS = 0.14;

% Degradation rate of 10.4%
degradation_sc = 10.4/100;
L_d_sc = (1 - degradation_sc).^(te/365.2422); % Lifetime degredation

%Power generated by solar array
gen_sc = eff_CIGS.*I.*area_sc.*L_d_sc;

% If considering
%{
% Can generate power between 90 -> 270 degrees, else turned off 
solar_on_sc = TA>=90 & TA <= 270 ;

% When operational
gen_operational_sc = gen_sc.*solar_on_sc.*abs(cosd(TA));
energy_gen_cum_sc = cumtrapz(te,gen_operational_sc);
%}

%
% When operational
gen_operational_sc = gen_sc.*abs(cosd(TA));
energy_gen_cum_sc = cumtrapz(te,gen_operational_sc);
%}

fig3 = figure(3);
hold on
plot(te/yr, energy_gen_cum_sc)
grid on
grid minor
box on
title('Cumulative Energy Generation - Solar Cells')
xlabel('Time elapsed [yrs]')
ylabel('Total energy generated [Wh]')
xlim([0 21])
hold off

fig3.Units = 'inches';
fig3.Position(3) = 6;
fig3.Position(4) = 3;
set(fig3.Children, 'FontName', 'Arial', 'FontSize', 11);
%print('SolarCellGen', '-depsc')

%% Energy balance

% total energy generation
tot_energy_gen = gen_operational_sc + gen_operational_sp;

% 2 kg of battery -> 600 Whrs of capacity stored 
battery_mass = 1.2; % in kg

capacity = zeros(te_maxtrack,1);
capacity(1) = battery_mass * LiCFx_energydensity;

% interpolate energy generation values to get values for each day
daily_tot_energy_gen = interp1q(te,tot_energy_gen,t_tracking);
daily_tot_energy_gen(6919:6929) = 0;

fig4 = figure(4);
hold on
plot(t_tracking/yr,daily_tot_energy_gen)
title('Energy Generation - Both sources')
grid on
grid minor
box on
xlabel('Time elapsed [yrs]')
ylabel('Energy generated [Wh]')
hold off

fig4.Units = 'inches';
fig4.Position(3) = 6;
fig4.Position(4) = 3;
set(fig4.Children, 'FontName', 'Arial', 'FontSize', 11);
%print('AllCumGen', '-depsc')

% Energy balance time!!! Wooohooooooooo!!! please make this end
for i = 2:te_maxtrack
    vari = capacity(i-1) + daily_tot_energy_gen(i) - energy_con(i,2);
    if vari >= capacity(1)
        capacity(i) = capacity(1);
    else
        capacity(i) = vari;
    end
    if capacity(i) <= 0
        %%print(['Out of energy! Happened on day',i,'.'])
        index_fail = i;
        break
    end
end

fig5 = figure(5);
hold on
plot(t_tracking/yr, capacity)
grid on
grid minor
box on
title('Energy capacity balance of sailcraft')
xlabel('Years elapsed [yrs]')
ylabel('Energy capacity in system [Wh]')
hold off

fig5.Units = 'inches';
fig5.Position(3) = 6;
fig5.Position(4) = 3;
set(fig5.Children, 'FontName', 'Arial', 'FontSize', 11);
%print('EnergyBalance', '-depsc')

battery_death_radius = interp1q(te,R,t_tracking(index_fail))/AU;

fig6 = figure(6);
hold on
plot(te/yr,R/AU,'-b')
plot(t_tracking(index_fail)/yr,battery_death_radius,'xr')
grid on
grid minor
box on
title('Sailcraft position over mission')
xlabel('Time elapsed [yrs]')
ylabel('Radius from Sun [AU]')
legend('Radius position of sailcraft','Death point','Location','Northwest')
xlim([0 21])
hold off

fig6.Units = 'inches';
fig6.Position(3) = 6;
fig6.Position(4) = 3;
set(fig6.Children, 'FontName', 'Arial', 'FontSize', 11);
%print('DeathPoint', '-depsc')

%% Temperature balance
rho_r = 0.91; % Reflectivitiy of reflective surface
rho_e = 0.16; % Reflectivity of emissive surface
SBc = 5.670374419e-8; % Steffan boltzmann constant
R_tracking = interp1q(te,R,t_tracking);
chi2_daily = interp1q(te,chi2,t_tracking);
chi3_daily = interp1q(te,chi3,t_tracking);
Temp_sail_daily = Temp_sail(R_tracking*1000,chi2_daily,chi3_daily,rho_r,rho_e);


fig7 = figure(7);
hold on
plot(t_tracking/yr,Temp_sail_daily,'-b')
grid on
grid minor
box on
title('Temperature of sailcraft over mission')
xlabel('Time elapsed [yrs]')
ylabel('Sailcraft temperature [K]')
xlim([0 21])
hold off

fig7.Units = 'inches';
fig7.Position(3) = 6;
fig7.Position(4) = 3;
set(fig7.Children, 'FontName', 'Arial', 'FontSize', 11);
%print('TempTime', '-depsc')