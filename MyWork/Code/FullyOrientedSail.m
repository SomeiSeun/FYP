% Full orientation control method can be utilised for complete analysis as
% previous analysis assumes perfect deployment such that if AoA is 90, then
% entire sail face is perfectly perpendicular to the Sun-Sail vector.
% Tanmay Ubgade | 220531

%% housekeeping
clear all
close all
clc

phi = pi/4; % Spin axis rotation
gam = rad2deg(45); % Width axis rotation (Ideally near 90 degrees)
psi = pi/4; % AoA

side = sqrt(28.6);

r_sail = ones(100).*linspace(-side/2,side/2); % x in body frame
w_sail = ones(100).*(linspace(side/2,-side/2)'); % y in body frame
h_sail = zeros(100);

R = RotZYZ(phi,gam,psi);
X = h_sail; Y = X; Z = X;
for i = 1:length(r_sail)
    for j = 1:length(r_sail)
        values = R*[r_sail(i,j);w_sail(i,j);h_sail(i,j)];
        X(i,j) = values(1);
        Y(i,j) = values(2); 
        Z(i,j) = values(3);
    end
end
