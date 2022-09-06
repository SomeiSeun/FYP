% Script to calculate the total power generation for sailcraft
% Tanmay Ubgade 220321

%% housekeeping
clear
clc

%% Subsystem weights in kg
TTC_weight = 2.2e-3;
CDH_weight = 0.1;

%% Power requirement
P_req = 5; % in watts
M_dry = 0.146; %

%% Solar Panel option

Area = 4*(0.28*0.08);

SPV_Weight = 0.04*P_req*4;

SPV_battery_weight = 35/45; 

remaining_weight = (0.02+0.025)*P_req + 0.02*M_dry ; % Power control unit, regulator/converter, wiring

% Total Mass
Mass_tot = SPV_battery_weight + SPV_Weight + remaining_weight