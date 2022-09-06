% Code to calculate the torque on the sail at a point in time
% By Tanmay Ubgade 220530

%% houskeeping
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
% Dimension of Sail
side = sqrt(28.6);

% Rotation of body about spin axis
%angle_spinaxis = 60; % Orientation about spin axis (clockwise)
%theta = 360-angle_spinaxis; % Anti-clockwise value
theta = 60;
% AoA of Sail to Sun, 0 degrees means parallel w Sun-Sail vector, 90
% degrees means perfectly facing the Sun
test = 45;
alpha = test; 

% Distance from Sun
[te, R, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ImportGMATData(1);
R_maxtrack = 8*AU;
maxtrack_index = find(R<R_maxtrack, 1, 'last');
te_maxtrack = floor(te(maxtrack_index));
t_tracking = 1:1:te_maxtrack; t_tracking = t_tracking';
R_tracking = interp1q(te,R,t_tracking);

radius = mean(R_tracking)*1000;
radius_mean_AU = radius/(1000*AU);
radius = AU*0.13*1000; % in meters


%% Moments of inertia calculated from Fusion 360

% Values taken about the origin in g/mm^2

% For Sail
Ixx = 3.335E+10;
Ixy = 1.138E-05;
Ixz = 1.349E-07;
Iyx = 1.138E-05;
Iyy = 3.335E+10;
Iyz = 1.349E-07;
Izx = 1.349E-07;
Izy = 1.349E-07;
Izz = 6.664E+10;

I_Sail = [Ixx, Ixy, Ixz; Iyx, Iyy, Iyz; Izx, Izy, Izz];
clear Ixx Ixy Ixz Iyx Iyy Iyz Izx Izy Izz

% For CubeSat
Ixx = 1.017E+08;
Ixy = 0.00;
Ixz = 0.00;
Iyx = 0.00;
Iyy = 1.017E+08;
Iyz = 0.00;
Izx = 0.00;
Izy = 0.00;
Izz = 5.500E+06;

I_CubeSat = [Ixx, Ixy, Ixz; Iyx, Iyy, Iyz; Izx, Izy, Izz];

I_tot = (I_Sail + I_CubeSat)*(1000^2)/1000;
 
%% SRP over sail(2D Case)

% Forming discretised nodes of x and y axis over the sail
x = ones(100).*linspace(-side/2,side/2); % X in body frame
y = ones(100).*(linspace(side/2,-side/2)'); % y in body frame

% Rotating coordinates based on orientation of Sail about spin axis
x_new = x.*cosd(theta)- y.*sind(theta);
y_new = y.*cosd(theta)+ x.*sind(theta);
% x and y axis is in same orientation as before, but nodes are mapped to
% their new position.


r = ones(100).*linspace(min([x_new,y_new],[],'all'),max([x_new,y_new],[],'all'))*cosd(alpha); % r in Sun-body frame

% Account for Angle of Attack

global_coords_r = radius + r;

SRP_local = SolarIntensity(global_coords_r).*sind(alpha)./c; % SRP across length of sail on each element
SRP_local_mu = SRP_local*10^6;

%
figure(3)
hold on 
surf(x_new, y_new, SRP_local_mu)
title('SRP over Sail - rotated frame')
xlabel('Radius axis [m]')
ylabel('Perp axis [m]')
zlabel('SRP [\muPa]')
grid on
view([6,6,2]) %(max(max(SRP_local_mu))+min(min(SRP_local_mu)))/2
hold off

figure(4)
hold on 
%plot3(x_new, y_new, SRP_local_mu, 'bo','MarkerSize',0.3)
surf(x, y, SRP_local_mu)
title('SRP over Sail - body frame')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('SRP [\muPa]')
grid on
view([6,6,2]) %(max(max(SRP_local_mu))+min(min(SRP_local_mu)))/2
hold off
%}

%% Element Sizes
elem_size = x(1,2) - x(1,1);
elem_area = elem_size^2;

%% Computing centroids of each element
% Need to consider 90 degree case due to way array is formed. else the
% resultant centroids and torques become 0. Converting angles within +-45
% degrees of 90 and 270 degrees helps account for these cases as properly
% travelling across matrix in various directions will prove to be quite
% complicated

x_centroids = zeros(length(x_new)-1);
y_centroids = x_centroids;
for i = 1:length(x_centroids) % rows
    for j = 1:length(x_centroids) % columns
        x_centroids(i,j) = mean(x_new(i:i+1,j:j+1),'all')/4;
        y_centroids(i,j) = mean(y_new(i:i+1,j:j+1),'all')/4;
    end
end

%% Computing SRP on each element
SRP_local_elem = zeros(length(x_centroids));
for i = 1:length(x_centroids) % rows
    for j = 1:length(x_centroids) % columns
        SRP_local_elem(i,j) = mean(SRP_local(i:i+1,j:j+1),'all')/4;
    end
end


% Computing force on each element
force_local = SRP_local_elem.*elem_area;
% Total force on sail accounting 
force_tot = sum(force_local,'all'); 

X_CEN_force_all = zeros(length(x_centroids),1);
Y_CEN_force_all = zeros(length(y_centroids),1);

%if (angle_spinaxis > 45 && angle_spinaxis < 135) || (angle_spinaxis > 225 && angle_spinaxis < 315)
%    for i = 1:length(x_centroids)
%        Y_CEN_force_all(i) = trapz(y_centroids(i,:),force_local(i,:))/force_tot;
%        X_CEN_force_all(i) = trapz(x_centroids(:,i),force_local(:,i))/force_tot;
%    end
%else
    for i = 1:length(x_centroids)
        X_CEN_force_all(i) = trapz(x_centroids(i,:),force_local(i,:))/force_tot;
        Y_CEN_force_all(i) = trapz(y_centroids(:,i),force_local(:,i))/force_tot;
    end    
%end

X_CEN_force = mean(X_CEN_force_all);
Y_CEN_force = mean(Y_CEN_force_all);
X_CG = 0.01;
Y_CG = 0.01;
%% Torque or Moment on Sail about each axis
T1 = force_tot*(X_CEN_force-X_CG);
T2 = force_tot*(Y_CEN_force-Y_CG);
% Axis 3 is spin axis where omega_dot_3 = 0 -> omega_3 = const.

%% Euler Equations to solve for angular velocity induced
I = diag(I_tot); % rest are negligible
I3 = I(3);
It = I(1);

roll_rate = sqrt((0.01*force_tot)/(1*It));

% As we are looking at initial induced angular velocity, omega1 and omega2
% are considered to be zero hence simplified equation.
omega_dot_1 = rad2deg(T1/It)
omega_dot_2 = rad2deg(T2/It)

omega_1_10d = omega_dot_1*(10*day*60*60)
omega_2_10d = omega_dot_2*(10*day*60*60)

angle_1_10d = (omega_dot_1*(10*day*60*60)^2)/2
angle_2_10d = (omega_dot_2*(10*day*60*60)^2)/2

toc