% Function that provides output for Torque on Sail
% Code by Tanmay Ubgade | 220530

function [T] = TorqueOnSailOffset(sideLength, spinAxRot, alpha, radius, offset)
%% Inputs
% sideLength - Side length of sail in m
% spinAxRot - Anticlockwise rotation about spin axis in degrees
% alpha - Angle of Attack of the Sail in degrees
% radius - Position relative to the Sun in m

%% Constants
% Speed of Light
c = 299792458; % [m/s]

%% SRP over sail(2D Case)

% Forming discretised nodes of x and y axis over the sail
x = ones(100).*linspace(-sideLength/2,sideLength/2); % x-axis in body frame
y = ones(100).*(linspace(sideLength/2,-sideLength/2)'); % y-axis in body frame

% Rotating coordinates based on orientation of Sail about spin axis
x_new = x.*cosd(spinAxRot) - y.*sind(spinAxRot);
y_new = y.*cosd(spinAxRot) + x.*sind(spinAxRot);
% x and y axis is in same orientation as before, but nodes are mapped to
% their new position.

% Radius axis to get variation of sail at radius point
r = ones(100).*linspace(min(x_new,[],'all'),max(x_new,[],'all'))*cosd(alpha); % r in Sun-body frame
% Convert variation of length to acquire Sun radius coordinates
global_coords_r = radius + r;

% Compute SRP over each element
SRP_local = SolarIntensity(global_coords_r)*sind(alpha)./c; % SRP across length of sail on each element

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
%% Force and force centroid computation
force_local = SRP_local_elem.*elem_area;
% Total force on sail accounting 
force_tot = sum(force_local,'all')*(1+0.6); % 0.6 is the reflectance factor 

X_CEN_force_all = zeros(length(x_centroids),1);
Y_CEN_force_all = zeros(length(y_centroids),1);

for i = 1:length(x_centroids)
    X_CEN_force_all(i) = trapz(x_centroids(i,:),force_local(i,:))/force_tot;
    Y_CEN_force_all(i) = trapz(y_centroids(:,i),force_local(:,i))/force_tot;
end

X_CEN_force = mean(X_CEN_force_all);
Y_CEN_force = mean(Y_CEN_force_all);
X_CG = offset;
Y_CG = offset;

%% Torque or Moment on Sail about each axis
T1 = force_tot*(X_CEN_force-X_CG);
T2 = force_tot*(Y_CEN_force-Y_CG);
T = [T1,T2,0]; % T = [T1, T2, T3], 1 is x, 2 is y, 3 is z (spin-axis)
end