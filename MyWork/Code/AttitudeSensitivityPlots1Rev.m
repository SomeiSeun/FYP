% Code providing plots for Attitude Sensitvity Script 

%% housekeeping
clear all
close all
clc

%% Import data
load('TorqueResult_1Rev.mat')

%% Inputs
row = 8; col = 4;

%% Plotting

PlotT1 = zeros(1,length(spinAxRot)); PlotT2 = PlotT1;
for i = 1:row
    % 1 is 45 deg, 2 is 135 deg, 3 is 225 deg, 4 is 315 deg
    for j = 1:col
        for k = 1:length(spinAxRot)
            PlotT1(k) =  T1(i,j,k);
            PlotT2(k) =  T2(i,j,k);
        end
        figure(1)
        hold on
        plot(spinAxRot, PlotT1);
        xlabel(['Spin Axis Rotation [' char(176) ']'])
        ylabel('Torque about x axis [Nm]')
        title('T_1 over one revolution')
        xlim([0,360])
        hold off
        % T2
        figure(2)
        hold on
        plot(spinAxRot, PlotT2);
        xlabel(['Spin Axis Rotation [' char(176) ']'])
        ylabel('Torque about y axis [Nm]')
        title('T_2 over one revolution')
        xlim([0,360])
        hold off
    end
end 

% Saving plots