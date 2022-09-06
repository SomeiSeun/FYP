% Code to verify that sail fits within CubeSat
% by Tanmay Ubgade | 220522

%% housekeeping
clear all
clc

%% Input parameters
areas = [7.5,12.2,28.3,454]; % areas in m^2

%% CubeSat variables in mm
length_CS = 400;
width_CS = 400;
height_CS = 300;
wall_thickness = 2;
corner_thickness = 5;
diag_length = sqrt(length_CS^2 + width_CS^2);

%% Folding values

sail_side = sqrt(areas); % Length of outmost part of sail in m
sail_edge = sqrt(2*sail_side.^2)*1000 - 10; % Edge of each triangle, converted to mm. Assuming 5mm spool centre

% 500 folds would be needed to make the sail thickness into 1mm thickness,
% let's assume 50 folds (excessive?) will be made, leading to 0.1mm
% thickness of the folded sail triangle
% This now is to be wrapped around the spool centre

%% Wrapping calculations
% Following calculations were made assuming that only one quadrant of the
% sail is being wrapped
%{
% Spool diam is 9 cm, available area for wrapping is 8 cm.
wrap_thickness = [0,0,0];
r = 5; % 5 mm, 1 cm diam
t = 0.2;
for i = 1:3
    l = sail_edge(i);
    length_rem = l;
    layers = 0;
    while length_rem > 0
        circ = 2*pi*(r + layers*t);
        length_rem = length_rem - circ;
        layers = layers + 1;
    end
    wrap_thickness(i) = 2*(layers*t + r);
end
empty_spool_space = 80 - wrap_thickness;
%}

% If wrapping all quadrants at the same time, each quadrant wraps a
% quarter-circumference before an additional layer is to be considered

wrap_thickness = [0,0,0,0];
r = 5; % 5 mm, 1 cm diam of the center of spool
t = 0.1; % Of one layer of sail that is to be wrapped
for i = 1:4
    l = sail_edge(i);
    length_rem = l;
    layers = 0;
    while length_rem > 0
        circ = 2*pi*(r + layers*t);
        length_rem = length_rem - circ/4;
        layers = layers + 1;
    end
    wrap_thickness(i) = 2*(layers*t + r);
end
% Spool outer diameter is distance between corners of CubeSat minus an
% additional 1 mm at each corner for clearance
spool_outer_diam = 384; 
empty_spool_space = spool_outer_diam - (wrap_thickness + 2*r);

%% Corner masses
mass_ball = 0.0789*[0.069, 0.087, 0.081244,1];
lead_density = 11.29/(1000*1000);
volume_ball = mass_ball/lead_density;
diam = 2*((0.75*volume_ball/pi).^(1/3));
empty_spool_space = empty_spool_space - diam 

%% Structural mass
topbot_vol = (wall_thickness+1)*length_CS*width_CS;
side_vol = wall_thickness*height_CS*width_CS;
support_vol = sqrt(5)*height_CS;
tot_vol = topbot_vol*2 + side_vol*4 + support_vol*4;
mass_structure = tot_vol*0.00000281 % vol * density in kg/mm^3