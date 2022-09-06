% Code that showcases the sailcraft remains spin-stabilised
% By Tanmay Ubgade | 220531

%% housekeeping
clear all
close all
clc

tic
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
side = sqrt(28.6); % Side length of sail in m
spinAxRot_initial = 0; % Rotation about spin-axis in degrees
alpha = 90; % Angle of Attack in degrees
I_tot = MomOfInertia(); % Moment of inertia of sail in kg/m^2

%% Acquiring positions of Sail when at alpha of 45 degrees
[te, R, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, G, ~, ~, ~, ~, ~, ~] = ImportGMATData(1);
TA = G(:,2); % True anomaly in degrees

R_maxtrack = 8*AU;
maxtrack_index = find(R<R_maxtrack, 1, 'last');
te_maxtrack = floor(te(maxtrack_index));
t_tracking = 1:1:te_maxtrack; t_tracking = t_tracking';
R_tracking = interp1q(te,R,t_tracking);
TA_tracking = interp1q(te,TA,t_tracking);

AoAofI = [45,135,225,315]; % Angles of Attack of interest
% Index of numbers before the stated angle, interpolate with index after
deg45  = [303,886,1443,2015,2640,3388,4444,6925];
deg135 = [451,1018,1565,2128,2743,2477,4541,6955];
deg225 = [91,667,1264,1863,2514,3287,4372,6896];
deg315 = [216,810,1385,1971,2609,3368,4432,6922];
deg0   = [260,853,1416,2002,2631,3380,4441,6926]; % Go backwards
deg180 = [2,556,1137,1711,2317,3010,3917,5695]; % Go backwards

R_ofI = zeros(length(deg45),6);
for i = 1:length(deg45)
    R_ofI(i,1) = interpolator(45,TA_tracking(deg45(i)),...
                                TA_tracking(deg45(i)+1),...
                                R_tracking(deg45(i)),...
                                R_tracking(deg45(i)+1)); % 45 degrees
    R_ofI(i,2) = interpolator(135,TA_tracking(deg135(i)),...
                                TA_tracking(deg135(i)+1),...
                                R_tracking(deg135(i)),...
                                R_tracking(deg135(i)+1)); % 135 degrees
    R_ofI(i,3) = interpolator(225,TA_tracking(deg225(i)),...
                                TA_tracking(deg225(i)+1),...
                                R_tracking(deg225(i)),...
                                R_tracking(deg225(i)+1)); % 225 degrees
    R_ofI(i,4) = interpolator(315,TA_tracking(deg315(i)),...
                                TA_tracking(deg315(i)+1),...
                                R_tracking(deg315(i)),...
                                R_tracking(deg315(i)+1)); % 315 degrees
    R_ofI(i,5) = interpolator(0,TA_tracking(deg0(i)-1),...
                                TA_tracking(deg0(i)),...
                                R_tracking(deg0(i)-1),...
                                R_tracking(deg0(i))); % 0 degrees
    R_ofI(i,6) = interpolator(0,TA_tracking(deg180(i)-1),...
                                TA_tracking(deg180(i)),...
                                R_tracking(deg180(i)-1),...
                                R_tracking(deg180(i))); % 0 degrees
end
R_ofI = R_ofI*1000; % Needs to be in meters
[row, col] = size(R_ofI);
%
%% Propagate attitude over 1 revolution
T1 = zeros(row, col); T2 = T1; %T3 = T1;

spinAxRot = 0:0.5:360;

% Choose which TorqueOnSail function to use depending on if it is offset
% case or idealised case

for i = 8%1:row
    % 1 is 45 deg, 2 is 135 deg, 3 is 225 deg, 4 is 315 deg
    for j = 5%5:col
        for k = 1:length(spinAxRot)
            Torque = TorqueOnSailIdeal(side, spinAxRot(k), alpha, R_ofI(i,j));
            T1(i,j,k) = Torque(1);
            T2(i,j,k) = Torque(2);
            %T3(i,j,k) = Torque(3);
            clear Torque
        end
    end
end
NetTorque1 = zeros(row, col); NetTorque2 = zeros(row, col);
%}
RPM = 0.63; t_rev = 60/RPM;
t_rot = linspace(0,t_rev,length(spinAxRot));
PlotT1 = zeros(1,length(spinAxRot)); PlotT2 = PlotT1;
for i = 8%1:row
    % 1 is 45 deg, 2 is 135 deg, 3 is 225 deg, 4 is 315 deg
    for j = 5%5:col
        for k = 1:length(spinAxRot)
            PlotT1(k) =  T1(i,j,k);
            PlotT2(k) =  T2(i,j,k);
        end
        NetTorque1(i,j) = trapz(t_rot,PlotT1);
        NetTorque2(i,j) = trapz(t_rot,PlotT2);
        fig1 = figure(1);
        hold on
        plot(spinAxRot, PlotT1);
        xlabel(['Spin Axis Rotation [',char(176),']'])
        ylabel('Torque about x axis [Nm]')
        title('T_1 over one revolution')
        xlim([0,360])
        grid minor
        box on
        yline(0)
        hold off
        % T2
        fig2 = figure(2);
        hold on
        plot(spinAxRot, PlotT2);
        xlabel(['Spin Axis Rotation [',char(176),']'])
        ylabel('Torque about y axis [Nm]')
        title('T_2 over one revolution')
        xlim([0,360])
        grid minor
        box on
        yline(0)
        hold off
    end
end

toc

fig1.Units = 'inches';
fig1.Position(3) = 3.5;
fig1.Position(4) = 2.8;
set(fig1.Children, 'FontName', 'Arial', 'FontSize', 11);
print('Torque_x_ideal', '-depsc')

fig2.Units = 'inches';
fig2.Position(3) = 3.5;
fig2.Position(4) = 2.8;
set(fig2.Children, 'FontName', 'Arial', 'FontSize', 11);
print('Torque_y_ideal', '-depsc')

%% Functions

function [y] = interpolator(x,x1,x2,y1,y2)
    y = y1 + (x-x1)*(y2-y1)/(x2-x1);
end
