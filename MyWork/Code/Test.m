clear
clc

side = 100; % 0 degree
diag = sqrt(2*side^2); % 45 degree

angle = 190;
theta = 360-angle;
dist = abs(side*cosd(theta) - side*sind(theta))
