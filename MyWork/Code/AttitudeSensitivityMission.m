% Code that simulates the accumulation and rotation of spin over the
% mission
% Tanmay Ubgade | 220601
%% housekeeping
clear all
close all
clc

%% Constants
% AU to km
AU = 149597870.7; % [km]
% Speed of Light
c = 299792458; % [m/s]
% Years to days
yrs = 365.2422; % [days]
% Hours in a day
day = 24;


%% Inputs
side = sqrt(28.6); % Side length of sail in m
I_tot = MomOfInertia(); % Moment of inertia of sail in kg/m^2
[te, R, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, G, ~, ~, ~, ~, ~, ~] = ImportGMATData(1);
TA = G(:,2); % True anomaly in ',char(176),'rees
te_s = te*day*60*60; % Time elapsed in seconds

% Spin-axis rotation
RPM = 0.63; 
t_rev = 60/RPM;
%% Torque Calculation loop
limit_21yrs = length(find((te < 21*yrs)));
T1 = zeros(1,limit_21yrs); T2 = T1; omega_dot_1 = T1; omega_dot_2 = T1;

offset = [-10, -5, -1, -0.1, -0.05, 0, 0.05, 0.1, 1, 5, 10]/1000; % CG offset value in meters

rot1 = zeros(1,length(offset)); rot2 = rot1; omega1 = rot1; omega2 = rot1;
count = 1;

for j = offset
    tic
    for i = 1:limit_21yrs
        spinAxRot = rem(te_s(i),t_rev)*360;
        Torque = TorqueOnSailOffset(side,spinAxRot,AoA(TA(i)),R(i),j);
        T1(i) = Torque(1);
        T2(i) = Torque(2);
        omega_dot = EulerEq12(Torque,I_tot);
        omega_dot_1(i) = omega_dot(1);
        omega_dot_2(i) = omega_dot(2);
        clear Torque omega_dot
    end

    %% Rotation rate
    omega1 = trapz(te_s(1:limit_21yrs), omega_dot_1);
    omega1_array = cumtrapz(te_s(1:limit_21yrs), omega_dot_1);

    omega2 = trapz(te_s(1:limit_21yrs), omega_dot_1);
    omega2_array = cumtrapz(te_s(1:limit_21yrs), omega_dot_2);

    %% Rotations
    rot1(count) = trapz(te_s(1:limit_21yrs), omega1_array);
    rot1_array = cumtrapz(te_s(1:limit_21yrs), omega1_array);

    rot2(count) = trapz(te_s(1:limit_21yrs), omega2_array);
    rot2_array = cumtrapz(te_s(1:limit_21yrs), omega2_array);

    %% Plots

    fig1 = figure(1);
    hold on 
    plot(te_s(1:limit_21yrs)/(day*yrs*60*60),rot1_array)
    title('X-axis rotation over mission')
    xlabel('Time elapsed [yrs]')
    ylabel(['X-axis rotation [',char(176),']'])
    grid on; grid minor
    box on
    xlim([0,21])
    hold off
    
    fig2 = figure(2);
    hold on 
    plot(te_s(1:limit_21yrs)/(day*yrs*60*60),rot2_array)
    title('Y-axis rotation over mission')
    xlabel('Time elapsed [yrs]')
    ylabel(['Y-axis rotation [',char(176),']'])
    grid on; grid minor
    box on
    xlim([0,21])
    hold off
    
    fig3 = figure(3);
    hold on
    plot(te_s(1:limit_21yrs)/(day*yrs*60*60),omega1_array)
    title('$\omega_1$ over mission','interpreter','latex')
    xlabel('Time elapsed [yrs]')
    ylabel(['X-axis angular velocity [',char(176),'s^{-1}]'])
    grid on; grid minor
    box on
    xlim([0,21])
    hold off

    fig4 = figure(4);
    hold on
    plot(te_s(1:limit_21yrs)/(day*yrs*60*60),omega2_array)
    title('$\omega_2$ over mission','interpreter','latex')
    xlabel('Time elapsed [yrs]')
    ylabel(['Y-axis angular velocity [',char(176),'s^{-1}]'])
    grid on; grid minor
    box on
    xlim([0,21])
    hold off
    
    fig5 = figure(5);
    hold on
    plot(te_s(1:limit_21yrs)/(day*yrs*60*60),omega_dot_1)
    title('$\dot{\omega}_1$ over mission','interpreter','latex')
    xlabel('Time elapsed [yrs]')
    ylabel(['X-axis angular acceleration [',char(176),'s^{-2}]'])
    grid on; grid minor
    box on
    xlim([0,21])
    hold off

    fig6 = figure(6);
    hold on 
    plot(te_s(1:limit_21yrs)/(day*yrs*60*60),omega_dot_2)
    title('$\dot{\omega}_2$ over mission','interpreter','latex')
    xlabel('Time elapsed [yrs]')
    ylabel(['Y-axis angular accleration [',char(176),'s^{-2}]'])
    grid on; grid minor
    box on
    xlim([0,21])
    hold off
    
    
    clear omega1_array omega2_array omega_dot_1 omega_dot_2 ...
          rot1_array rot2_array
    count = count + 1;
    toc
end


    
fig1.Units = 'inches';
fig1.Position(3) = 6;
fig1.Position(4) = 3;
set(fig1.Children, 'FontName', 'Arial', 'FontSize', 11);

fig2.Units = 'inches';
fig2.Position(3) = 6;
fig2.Position(4) = 3;
set(fig2.Children, 'FontName', 'Arial', 'FontSize', 11);

fig3.Units = 'inches';
fig3.Position(3) = 6;
fig3.Position(4) = 3;
set(fig3.Children, 'FontName', 'Arial', 'FontSize', 11);

fig4.Units = 'inches';
fig4.Position(3) = 6;
fig4.Position(4) = 3;
set(fig4.Children, 'FontName', 'Arial', 'FontSize', 11);

fig5.Units = 'inches';
fig5.Position(3) = 6;
fig5.Position(4) = 3;
set(fig5.Children, 'FontName', 'Arial', 'FontSize', 11);



fig6.Units = 'inches';
fig6.Position(3) = 6;
fig6.Position(4) = 3;
set(fig6.Children, 'FontName', 'Arial', 'FontSize', 11);



figure(1)
legend('y_{cg}= -10mm','y_{cg}= -5mm','y_{cg}= -1mm','y_{cg}= -0.1mm',...
        'y_{cg}= -0.05mm','y_{cg}= 0mm','y_{cg}= 0.05mm','y_{cg}= 0.1mm',...
        'y_{cg}= 1mm','y_{cg}= 5mm','y_{cg}= 10mm',...
        'Location','eastoutside','FontSize',8)
print('rot1', '-depsc')
figure(2)
legend('x_{cg}= -10mm','x_{cg}= -5mm','x_{cg}= -1mm','x_{cg}= -0.1mm',...
        'x_{cg}= -0.05mm','x_{cg}= 0mm','x_{cg}= 0.05mm','x_{cg}= 0.1mm',...
        'x_{cg}= 1mm','x_{cg}= 5mm','x_{cg}= 10mm',...
        'Location','eastoutside','FontSize',8)
print('rot2', '-depsc')
figure(3)
legend('y_{cg}= -10mm','y_{cg}= -5mm','y_{cg}= -1mm','y_{cg}= -0.1mm',...
        'y_{cg}= -0.05mm','y_{cg}= 0mm','y_{cg}= 0.05mm','y_{cg}= 0.1mm',...
        'y_{cg}= 1mm','y_{cg}= 5mm','y_{cg}= 10mm',...
        'Location','eastoutside','FontSize',8)
print('omega1', '-depsc')
figure(4)
legend('x_{cg}= -10mm','x_{cg}= -5mm','x_{cg}= -1mm','x_{cg}= -0.1mm',...
        'x_{cg}= -0.05mm','x_{cg}= 0mm','x_{cg}= 0.05mm','x_{cg}= 0.1mm',...
        'x_{cg}= 1mm','x_{cg}= 5mm','x_{cg}= 10mm',...
        'Location','eastoutside','FontSize',8)
print('omega2', '-depsc')
figure(5)
legend('y_{cg}= -10mm','y_{cg}= -5mm','y_{cg}= -1mm','y_{cg}= -0.1mm',...
        'y_{cg}= -0.05mm','y_{cg}= 0mm','y_{cg}= 0.05mm','y_{cg}= 0.1mm',...
        'y_{cg}= 1mm','y_{cg}= 5mm','y_{cg}= 10mm',...
        'Location','eastoutside','FontSize',8)
print('omegaDot1', '-depsc')
figure(6)
legend('x_{cg}= -10mm','x_{cg}= -5mm','x_{cg}= -1mm','x_{cg}= -0.1mm',...
        'x_{cg}= -0.05mm','x_{cg}= 0mm','x_{cg}= 0.05mm','x_{cg}= 0.1mm',...
        'x_{cg}= 1mm','x_{cg}= 5mm','x_{cg}= 10mm',...
        'Location','eastoutside','FontSize',8)
print('omegaDot2', '-depsc')

%% Functions
function alpha = AoA(TA)
    if (TA > 0) & (TA < 90)
        alpha = 90 - TA;
    elseif (TA > 90)  & (TA < 180)
        alpha = TA - 90;
    elseif (TA > 180) & (TA < 270)
        alpha = 180 - TA;
    else
        alpha = TA - 270;
    end
end