% This code imports, tabulates and plots data from GMAT orbital simulations
% for the trajectory of 'Sunbeam' to determine the critical perihelion
% for sun escape
% Written by Filippos Geragidis on 20/3/2021

%% Housekeeping
clear;
clc;
close all;


% I/O & preliminary processing

% Conversion factors and physical properties

% AU to km
AU = 149597870.7; % [km]
% Years to days
yr = 365.2422; % [days]
% Solar constants
I0 = 1365.4; % Solar intensity [W/m^2]
Cr = [0.91,0.16]; % Reflectivity coefficient vector for the two
% Solar gravitational constant [m^3/s^2]
mu = 1.32712440018e20;
% sail surfaces: [reflective,emissive] [-]
SBc = 5.670374419e-8; % Stefan-Boltzmann constant [W/(m^2 K^4)]
% Perihelion radius set for simulations [AU]
mass = [2.27,3,5,10,1.5]; % initial yaw angle


% Import data

% Latest simulation number
CurrentSimNo = 1;

% Storage counter variable - independent from loop variable
% in case of skipped datasets
n = 1;

% Import loop for each output dataset
for i = 1:CurrentSimNo
    
    % Generate filenames for each simulation
    a = ['Sim',num2str(i),'a.txt']; % Simxa.txt contains data on elapsed time, radial distance and SMA for sim number x
    b = ['Sim',num2str(i),'b.txt']; % Simxb.txt contains data on velocity, and 3-axis acceleration for sim number x
    c = ['Sim',num2str(i),'c.txt']; % Simxc.txt contains data on Euler angles for sim number x
    d = ['Sim',num2str(i),'d.txt']; % Simxd.txt contains data on Euler angle rates for sim number x
    e = ['Sim',num2str(i),'e.txt']; % Simxe.txt contains data on X, Y, and Z coordinates for sim number x
    
    % Import data from files
    
    % Open documents
    aID = fopen(a,'r');
    bID = fopen(b,'r');
    cID = fopen(c,'r');
    dID = fopen(d,'r');
    eID = fopen(e,'r');
    
    % Skip headers
    fscanf(aID,'%s %s %s %s',[4,1]);
    fscanf(bID,'%s %s %s %s',[4,1]);
    fscanf(cID,'%s %s %s',[3,1]);
    fscanf(dID,'%s %s %s',[3,1]);
    fscanf(eID,'%s %s %s',[3,1]);
    
    % Store data into file-specific matrices
    A = transpose(fscanf(aID,'%f %*f %f %f',[3,inf]));
    B = transpose(fscanf(bID,'%f %f %f %f',[4,inf]));
    C = transpose(fscanf(cID,' %f %f %f',[3,inf]));
    D = transpose(fscanf(dID,' %f %f %f',[3,inf]));
    E = transpose(fscanf(eID,' %f %f %f',[3,inf]));
    
    % Close files
    fclose(aID);
    fclose(bID);
    fclose(cID);
    fclose(dID);
    fclose(eID);
    
    % Store results into sim-specific vectors to avoid data size mismatches
    % Note that trajectories are separated into extrasolar and non,
    % in order to discard extraneous data past the maximum separation
    % point
    if n == 1
        if max(A(:,2)/AU) >= 122
            te1(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),1)/yr; % elapsed time [yrs] (converted from days)
            R1(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),2)/AU; % orbital distance [AU] (converted from km)
            SMA1(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),3)/AU; % Semi-major axis [AU] (converted from km)
            V1(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),1); % Velocity [km/s]
            Ax1(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),2)*1e3; % Normal acceleration [m/s]
            Ay1(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),3)*1e3; % Lateral acceleration [m/s]
            Az1(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),4)*1e3; % Vertical acceleration [m/s]
            E11(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),1); % Euler angle 1 [deg]
            E21(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),2); % Euler angle 2 [deg]
            E31(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),3); % Euler angle 3 [deg]
            E1r1(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),1); % Euler angle rate 1 [deg]
            E2r1(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),2); % Euler angle rate 2 [deg]
            E3r1(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),3); % Euler angle rate 3 [deg]
            X1(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),1); % X coordinate in Sun-fixed inertial coordinates [km]
            Y1(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),2); % Y coordinate in Sun-fixed inertial coordinates [km]
            Z1(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),3); % Z coordinate in Sun-fixed inertial coordinates [km]
            % Calculate sun-spacecraft angles and instantaneous spacecraft
            % temperatures
            for m = 1:length(X1) % for each data point
                chi21(m,1) = atand(Z1(m,1)/X1(m,1)) + E21(m,1); % inertial plane angle to body lateral plane
                chi31(m,1) = atand(Y1(m,1)/X1(m,1))+E31(m,1); % inertial lateral plane angle to body lateral plane
                % determining whether the sun-facing side is reflective
                % (front) or emmissive (rear)
                if chi21(m,1) < 90 && chi21(m,1) > -90 % front
                    if chi31(m,1) < 90 && chi31(m,1) > -90 % front
                        Cri = 1;
                    else % rear
                        Cri = 2;
                    end
                else % rear
                    if chi31(m,1) < 90 && chi31(m,1) > -90 % front
                        Cri = 2;
                    else % front
                        Cri = 1;
                    end
                end
                % emmissive equilibrium temperature
                T1(m,1) = (((1-Cr(Cri))*I0*abs(sind(chi21(m,1))*sind(chi31(m,1))))/((2-(Cr(1)+Cr(2)))*SBc*R1(m,1)))^0.25;
            end
            
            % Overall orbit properties
            Tmax(1,n) = max(T1); % Maximum instantaneous temperature [K]
            [qq, pp] = max(A(:,2));
            th(1,n) = A(pp,1); % Time to reach heliopause [yrs]
            Ve(1,n) = V1(end,1); % Exit velocity [km/s]
            Amax(1,n) = max([max(abs(Ax1)),max(abs(Ay1)),max(abs(Az1))]); % Max acceleration [m/s^2]
            Rmin(1,n) = min(R1); % Minimum solar flyby distance [AU]
        else
            [q, p1] = max(A(:,2));
            te1(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),1)/yr; % elapsed time [yrs] (converted from days)
            R1(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),2)/AU; % orbital distance [AU] (converted from km)
            SMA1(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),3)/AU; % Semi-major axis [AU] (converted from km)
            V1(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),1); % Velocity [km/s]
            Ax1(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),2)*1e3; % Normal acceleration [m/s]
            Ay1(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),3)*1e3; % Lateral acceleration [m/s]
            Az1(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),4)*1e3; % Vertical acceleration [m/s]
            E11(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),1); % Euler angle 1 [deg]
            E21(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),2); % Euler angle 2 [deg]
            E31(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),3); % Euler angle 3 [deg]
            E1r1(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),1); % Euler angle rate 1 [deg]
            E2r1(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),2); % Euler angle rate 2 [deg]
            E3r1(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),3); % Euler angle rate 3 [deg]
            X1(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),1); % X coordinate in Sun-fixed inertial coordinates [km]
            Y1(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),2); % Y coordinate in Sun-fixed inertial coordinates [km]
            Z1(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),3); % Z coordinate in Sun-fixed inertial coordinates [km]
            % Calculate sun-spacecraft angles and instantaneous spacecraft
            % temperatures
            for m = 1:length(X1) % for each data point
                chi21(m,1) = atand(Z1(m,1)/X1(m,1)) + E21(m,1); % inertial plane angle to body lateral plane
                chi31(m,1) = atand(Y1(m,1)/X1(m,1))+E31(m,1); % inertial lateral plane angle to body lateral plane
                % determining whether the sun-facing side is reflective
                % (front) or emmissive (rear)
                if chi21(m,1) < 90 && chi21(m,1) > -90 % front
                    if chi31(m,1) < 90 && chi31(m,1) > -90 % front
                        Cri = 1;
                    else % rear
                        Cri = 2;
                    end
                else % rear
                    if chi31(m,1) < 90 && chi31(m,1) > -90 % front
                        Cri = 2; % rear
                    else % front
                        Cri = 1; % front
                    end
                end
                % emmissive equilibrium temperature
                T1(m,1) = (((1-Cr(Cri))*I0*abs(sind(chi21(m,1))*sind(chi31(m,1))))/((2-(Cr(1)+Cr(2)))*SBc*R1(m,1)))^0.25;
            end
            
            % Overall orbit properties
            Tmax(1,n) = max(T1); % Maximum instantaneous temperature [K]
            th(1,n) = NaN; % Time to reach heliopause [yrs]
            Ve(1,n) = NaN; % Exit velocity [km/s]
            Amax(1,n) = max([max(abs(Ax1)),max(abs(Ay1)),max(abs(Az1))]); % Max acceleration [m/s^2]
            Rmin(1,n) = min(R1); % Minimum solar flyby distance [AU]
        end
    elseif n == 2
        if max(A(:,2)/AU) >= 120
            te2(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),1)/yr; % elapsed time [yrs] (converted from days)
            R2(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),2)/AU; % orbital distance [AU] (converted from km)
            SMA2(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),3)/AU; % Semi-major axis [AU] (converted from km)
            V2(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),1); % Velocity [km/s]
            Ax2(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),2)*1e3; % Normal acceleration [m/s]
            Ay2(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),3)*1e3; % Lateral acceleration [m/s]
            Az2(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),4)*1e3; % Vertical acceleration [m/s]
            E12(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),1); % Euler angle 1 [deg]
            E22(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),2); % Euler angle 2 [deg]
            E32(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),3); % Euler angle 3 [deg]
            E1r2(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),1); % Euler angle rate 1 [deg]
            E2r2(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),2); % Euler angle rate 2 [deg]
            E3r2(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),3); % Euler angle rate 3 [deg]
            X2(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),1); % X coordinate in Sun-fixed inertial coordinates [km]
            Y2(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),2); % Y coordinate in Sun-fixed inertial coordinates [km]
            Z2(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),3); % Z coordinate in Sun-fixed inertial coordinates [km]
            % Calculate sun-spacecraft angles and instantaneous spacecraft
            % temperatures
            for m = 1:length(X2) % for each data point
                chi22(m,1) = atand(Z2(m,1)/X2(m,1)) + E22(m,1); % inertial plane angle to body lateral plane
                chi32(m,1) = atand(Y2(m,1)/X2(m,1))+E32(m,1); % inertial lateral plane angle to body lateral plane
                % determining whether the sun-facing side is reflective
                % (front) or emmissive (rear)
                if chi22(m,1) < 90 && chi22(m,1) > -90 % front
                    if chi32(m,1) < 90 && chi32(m,1) > -90 % front
                        Cri = 1;
                    else % rear
                        Cri = 2;
                    end
                else % rear
                    if chi32(m,1) < 90 && chi32(m,1) > -90 % front
                        Cri = 2;
                    else % front
                        Cri = 1;
                    end
                end
                % emmissive equilibrium temperature
                T2(m,1) = (((1-Cr(Cri))*I0*abs(sind(chi22(m,1))*sind(chi32(m,1))))/((2-(Cr(1)+Cr(2)))*SBc*R2(m,1)))^0.25;
            end
            
            % Overall orbit properties
            Tmax(1,n) = max(T2); % Maximum instantaneous temperature [K]
            [qq, pp] = max(A(:,2));
            th(1,n) = A(pp,1); % Time to reach heliopause [yrs]
            Ve(1,n) = V2(end,1); % Exit velocity [km/s]
            Amax(1,n) = max([max(abs(Ax2)),max(abs(Ay2)),max(abs(Az2))]); % Max acceleration [m/s^2]
            Rmin(1,n) = min(R2); % Minimum solar flyby distance [AU]
        else
            [q, p2] = max(A(:,2));
            te2(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),1)/yr; % elapsed time [yrs] (converted from days)
            R2(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),2)/AU; % orbital distance [AU] (converted from km)
            SMA2(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),3)/AU; % Semi-major axis [AU] (converted from km)
            V2(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),1); % Velocity [km/s]
            Ax2(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),2)*1e3; % Normal acceleration [m/s]
            Ay2(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),3)*1e3; % Lateral acceleration [m/s]
            Az2(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),4)*1e3; % Vertical acceleration [m/s]
            E12(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),1); % Euler angle 1 [deg]
            E22(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),2); % Euler angle 2 [deg]
            E32(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),3); % Euler angle 3 [deg]
            E1r2(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),1); % Euler angle rate 1 [deg]
            E2r2(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),2); % Euler angle rate 2 [deg]
            E3r2(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),3); % Euler angle rate 3 [deg]
            X2(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),1); % X coordinate in Sun-fixed inertial coordinates [km]
            Y2(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),2); % Y coordinate in Sun-fixed inertial coordinates [km]
            Z2(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),3); % Z coordinate in Sun-fixed inertial coordinates [km]
            % Calculate sun-spacecraft angles and instantaneous spacecraft
            % temperatures
            for m = 1:length(X2) % for each data point
                chi22(m,1) = atand(Z2(m,1)/X2(m,1)) + E22(m,1); % inertial plane angle to body lateral plane
                chi32(m,1) = atand(Y2(m,1)/X2(m,1))+E32(m,1); % inertial lateral plane angle to body lateral plane
                % determining whether the sun-facing side is reflective
                % (front) or emmissive (rear)
                if chi22(m,1) < 90 && chi22(m,1) > -90 % front
                    if chi32(m,1) < 90 && chi32(m,1) > -90 % front
                        Cri = 1;
                    else % rear
                        Cri = 2;
                    end
                else % rear
                    if chi32(m,1) < 90 && chi32(m,1) > -90 % front
                        Cri = 2;
                    else % front
                        Cri = 1;
                    end
                end
                % emmissive equilibrium temperature
                T2(m,1) = (((1-Cr(Cri))*I0*abs(sind(chi22(m,1))*sind(chi32(m,1))))/((2-(Cr(1)+Cr(2)))*SBc*R2(m,1)))^0.25;
            end
            
            % Overall orbit properties
            Tmax(1,n) = max(T2); % Maximum instantaneous temperature [K]
            th(1,n) = NaN; % Time to reach heliopause [yrs]
            Ve(1,n) = NaN; % Exit velocity [km/s]
            Amax(1,n) = max([max(abs(Ax2)),max(abs(Ay2)),max(abs(Az2))]); % Max acceleration [m/s^2]
            Rmin(1,n) = min(R2); % Minimum solar flyby distance [AU]
        end
        
    elseif n == 3
        if max(A(:,2)/AU) >= 120
            te3(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),1)/yr; % elapsed time [yrs] (converted from days)
            R3(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),2)/AU; % orbital distance [AU] (converted from km)
            SMA3(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),3)/AU; % Semi-major axis [AU] (converted from km)
            V3(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),1); % Velocity [km/s]
            Ax3(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),2)*1e3; % Normal acceleration [m/s]
            Ay3(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),3)*1e3; % Lateral acceleration [m/s]
            Az3(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),4)*1e3; % Vertical acceleration [m/s]
            E13(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),1); % Euler angle 1 [deg]
            E23(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),2); % Euler angle 2 [deg]
            E33(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),3); % Euler angle 3 [deg]
            E1r3(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),1); % Euler angle rate 1 [deg]
            E2r3(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),2); % Euler angle rate 2 [deg]
            E3r3(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),3); % Euler angle rate 3 [deg]
            X3(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),1); % X coordinate in Sun-fixed inertial coordinates [km]
            Y3(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),2); % Y coordinate in Sun-fixed inertial coordinates [km]
            Z3(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),3); % Z coordinate in Sun-fixed inertial coordinates [km]
            % Calculate sun-spacecraft angles and instantaneous spacecraft
            % temperatures
            for m = 1:length(X3) % for each data point
                chi23(m,1) = atand(Z3(m,1)/X3(m,1)) + E23(m,1); % inertial plane angle to body lateral plane
                chi33(m,1) = atand(Y3(m,1)/X3(m,1))+E33(m,1); % inertial lateral plane angle to body lateral plane
                % determining whether the sun-facing side is reflective
                % (front) or emmissive (rear)
                if chi23(m,1) < 90 && chi23(m,1) > -90 % front
                    if chi33(m,1) < 90 && chi33(m,1) > -90 % front
                        Cri = 1;
                    else % rear
                        Cri = 2;
                    end
                else % rear
                    if chi33(m,1) < 90 && chi33(m,1) > -90 % front
                        Cri = 2;
                    else % front
                        Cri = 1;
                    end
                end
                % emmissive equilibrium temperature
                T3(m,1) = (((1-Cr(Cri))*I0*abs(sind(chi23(m,1))*sind(chi33(m,1))))/((2-(Cr(1)+Cr(2)))*SBc*R3(m,1)))^0.25;
            end
            
            % Overall orbit properties
            Tmax(1,n) = max(T3); % Maximum instantaneous temperature [K]
            [qq, pp] = max(A(:,2));
            th(1,n) = A(pp,1); % Time to reach heliopause [yrs]
            Ve(1,n) = V3(end,1); % Exit velocity [km/s]
            Amax(1,n) = max([max(abs(Ax3)),max(abs(Ay3)),max(abs(Az3))]); % Max acceleration [m/s^2]
            Rmin(1,n) = min(R3); % Minimum solar flyby distance [AU]
        else
            [q, p3] = max(A(:,2));
            te3(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),1)/yr; % elapsed time [yrs] (converted from days)
            R3(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),2)/AU; % orbital distance [AU] (converted from km)
            SMA3(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),3)/AU; % Semi-major axis [AU] (converted from km)
            V3(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),1); % Velocity [km/s]
            Ax3(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),2)*1e3; % Normal acceleration [m/s]
            Ay3(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),3)*1e3; % Lateral acceleration [m/s]
            Az3(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),4)*1e3; % Vertical acceleration [m/s]
            E13(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),1); % Euler angle 1 [deg]
            E23(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),2); % Euler angle 2 [deg]
            E33(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),3); % Euler angle 3 [deg]
            E1r3(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),1); % Euler angle rate 1 [deg]
            E2r3(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),2); % Euler angle rate 2 [deg]
            E3r3(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),3); % Euler angle rate 3 [deg]
            X3(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),1); % X coordinate in Sun-fixed inertial coordinates [km]
            Y3(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),2); % Y coordinate in Sun-fixed inertial coordinates [km]
            Z3(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),3); % Z coordinate in Sun-fixed inertial coordinates [km]
            % Calculate sun-spacecraft angles and instantaneous spacecraft
            % temperatures
            for m = 1:length(X3) % for each data point
                chi23(m,1) = atand(Z3(m,1)/X3(m,1)) + E23(m,1); % inertial plane angle to body lateral plane
                chi33(m,1) = atand(Y3(m,1)/X3(m,1))+E33(m,1); % inertial lateral plane angle to body lateral plane
                % determining whether the sun-facing side is reflective
                % (front) or emmissive (rear)
                if chi23(m,1) < 90 && chi23(m,1) > -90 % front
                    if chi33(m,1) < 90 && chi33(m,1) > -90 % front
                        Cri = 1;
                    else % rear
                        Cri = 2;
                    end
                else % rear
                    if chi33(m,1) < 90 && chi33(m,1) > -90 % front
                        Cri = 2;
                    else % front
                        Cri = 1;
                    end
                end
                % emmissive equilibrium temperature
                T3(m,1) = (((1-Cr(Cri))*I0*abs(sind(chi23(m,1))*sind(chi33(m,1))))/((2-(Cr(1)+Cr(2)))*SBc*R3(m,1)))^0.25;
            end
            
            % Overall orbit properties
            Tmax(1,n) = max(T3); % Maximum instantaneous temperature [K]
            th(1,n) = NaN; % Time to reach heliopause [yrs]
            Ve(1,n) = NaN; % Exit velocity [km/s]
            Amax(1,n) = max([max(abs(Ax3)),max(abs(Ay3)),max(abs(Az3))]); % Max acceleration [m/s^2]
            Rmin(1,n) = min(R3); % Minimum solar flyby distance [AU]
        end
    elseif n == 4
        if max(A(:,2)/AU) >= 120
            te4(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),1)/yr; % elapsed time [yrs] (converted from days)
            R4(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),2)/AU; % orbital distance [AU] (converted from km)
            SMA4(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),3)/AU; % Semi-major axis [AU] (converted from km)
            V4(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),1); % Velocity [km/s]
            Ax4(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),2)*1e3; % Normal acceleration [m/s]
            Ay4(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),3)*1e3; % Lateral acceleration [m/s]
            Az4(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),4)*1e3; % Vertical acceleration [m/s]
            E14(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),1); % Euler angle 1 [deg]
            E24(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),2); % Euler angle 2 [deg]
            E34(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),3); % Euler angle 3 [deg]
            E1r4(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),1); % Euler angle rate 1 [deg]
            E2r4(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),2); % Euler angle rate 2 [deg]
            E3r4(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),3); % Euler angle rate 3 [deg]
            X4(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),1); % X coordinate in Sun-fixed inertial coordinates [km]
            Y4(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),2); % Y coordinate in Sun-fixed inertial coordinates [km]
            Z4(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),3); % Z coordinate in Sun-fixed inertial coordinates [km]
            % Calculate sun-spacecraft angles and instantaneous spacecraft
            % temperatures
            for m = 1:length(X4) % for each data point
                chi24(m,1) = atand(Z4(m,1)/X4(m,1)) + E24(m,1); % inertial plane angle to body lateral plane
                chi34(m,1) = atand(Y4(m,1)/X4(m,1))+E34(m,1); % inertial lateral plane angle to body lateral plane
                % determining whether the sun-facing side is reflective
                % (front) or emmissive (rear)
                if chi24(m,1) < 90 && chi24(m,1) > -90 % front
                    if chi34(m,1) < 90 && chi34(m,1) > -90 % front
                        Cri = 1;
                    else % rear
                        Cri = 2;
                    end
                else % rear
                    if chi34(m,1) < 90 && chi34(m,1) > -90 % front
                        Cri = 2;
                    else % front
                        Cri = 1;
                    end
                end
                % emmissive equilibrium temperature
                T4(m,1) = (((1-Cr(Cri))*I0*abs(sind(chi24(m,1))*sind(chi34(m,1))))/((2-(Cr(1)+Cr(2)))*SBc*R4(m,1)))^0.25;
            end
            
            % Overall orbit properties
            Tmax(1,n) = max(T4); % Maximum instantaneous temperature [K]
            [qq, pp] = max(A(:,2));
            th(1,n) = A(pp,1); % Time to reach heliopause [yrs]
            Ve(1,n) = V4(end,1); % Exit velocity [km/s]
            Amax(1,n) = max([max(abs(Ax4)),max(abs(Ay4)),max(abs(Az4))]); % Max acceleration [m/s^2]
            Rmin(1,n) = min(R4); % Minimum solar flyby distance [AU]
        else
            [q, p4] = max(A(:,2));
            te4(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),1)/yr; % elapsed time [yrs] (converted from days)
            R4(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),2)/AU; % orbital distance [AU] (converted from km)
            SMA4(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),3)/AU; % Semi-major axis [AU] (converted from km)
            V4(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),1); % Velocity [km/s]
            Ax4(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),2)*1e3; % Normal acceleration [m/s]
            Ay4(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),3)*1e3; % Lateral acceleration [m/s]
            Az4(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),4)*1e3; % Vertical acceleration [m/s]
            E14(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),1); % Euler angle 1 [deg]
            E24(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),2); % Euler angle 2 [deg]
            E34(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),3); % Euler angle 3 [deg]
            E1r4(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),1); % Euler angle rate 1 [deg]
            E2r4(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),2); % Euler angle rate 2 [deg]
            E3r4(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),3); % Euler angle rate 3 [deg]
            X4(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),1); % X coordinate in Sun-fixed inertial coordinates [km]
            Y4(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),2); % Y coordinate in Sun-fixed inertial coordinates [km]
            Z4(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),3); % Z coordinate in Sun-fixed inertial coordinates [km]
            % Calculate sun-spacecraft angles and instantaneous spacecraft
            % temperatures
            for m = 1:length(X4) % for each data point
                chi24(m,1) = atand(Z4(m,1)/X4(m,1)) + E24(m,1); % inertial plane angle to body lateral plane
                chi34(m,1) = atand(Y4(m,1)/X4(m,1))+E34(m,1); % inertial lateral plane angle to body lateral plane
                % determining whether the sun-facing side is reflective
                % (front) or emmissive (rear)
                if chi24(m,1) < 90 && chi24(m,1) > -90 % front
                    if chi34(m,1) < 90 && chi34(m,1) > -90 % front
                        Cri = 1;
                    else % rear
                        Cri = 2;
                    end
                else % rear
                    if chi34(m,1) < 90 && chi34(m,1) > -90 % front
                        Cri = 2;
                    else % front
                        Cri = 1;
                    end
                end
                % emmissive equilibrium temperature
                T4(m,1) = (((1-Cr(Cri))*I0*abs(sind(chi24(m,1))*sind(chi34(m,1))))/((2-(Cr(1)+Cr(2)))*SBc*R4(m,1)))^0.25;
            end
        end
            
            % Overall orbit properties
            Tmax(1,n) = max(T4); % Maximum instantaneous temperature [K]
            th(1,n) = NaN; % Time to reach heliopause [yrs]
            Ve(1,n) = NaN; % Exit velocity [km/s]
            Amax(1,n) = max([max(abs(Ax4)),max(abs(Ay4)),max(abs(Az4))]); % Max acceleration [m/s^2]
            Rmin(1,n) = min(R4); % Minimum solar flyby distance [AU]

    elseif n == 5
        if max(A(:,2)/AU) >= 120
            te5(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),1)/yr; % elapsed time [yrs] (converted from days)
            R5(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),2)/AU; % orbital distance [AU] (converted from km)
            SMA5(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),3)/AU; % Semi-major axis [AU] (converted from km)
            V5(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),1); % Velocity [km/s]
            Ax5(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),2)*1e3; % Normal acceleration [m/s]
            Ay5(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),3)*1e3; % Lateral acceleration [m/s]
            Az5(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),4)*1e3; % Vertical acceleration [m/s]
            E15(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),1); % Euler angle 1 [deg]
            E25(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),2); % Euler angle 2 [deg]
            E35(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),3); % Euler angle 3 [deg]
            E1r5(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),1); % Euler angle rate 1 [deg]
            E2r5(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),2); % Euler angle rate 2 [deg]
            E3r5(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),3); % Euler angle rate 3 [deg]
            X5(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),1); % X coordinate in Sun-fixed inertial coordinates [km]
            Y5(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),2); % Y coordinate in Sun-fixed inertial coordinates [km]
            Z5(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),3); % Z coordinate in Sun-fixed inertial coordinates [km]
            % Calculate sun-spacecraft angles and instantaneous spacecraft
            % temperatures
            for m = 1:length(X5) % for each data point
                chi25(m,1) = atand(Z5(m,1)/X5(m,1)) + E25(m,1); % inertial plane angle to body lateral plane
                chi35(m,1) = atand(Y5(m,1)/X5(m,1))+E35(m,1); % inertial lateral plane angle to body lateral plane
                % determining whether the sun-facing side is reflective
                % (front) or emmissive (rear)
                if chi25(m,1) < 90 && chi25(m,1) > -90 % front
                    if chi35(m,1) < 90 && chi35(m,1) > -90 % front
                        Cri = 1;
                    else % rear
                        Cri = 2;
                    end
                else % rear
                    if chi35(m,1) < 90 && chi35(m,1) > -90 % front
                        Cri = 2;
                    else % front
                        Cri = 1;
                    end
                end
                % emmissive equilibrium temperature
                T5(m,1) = (((1-Cr(Cri))*I0*abs(sind(chi25(m,1))*sind(chi35(m,1))))/((2-(Cr(1)+Cr(2)))*SBc*R5(m,1)))^0.25;
            end
            
            % Overall orbit properties
            Tmax(1,n) = max(T5); % Maximum instantaneous temperature [K]
            [qq, pp] = max(A(:,2));
            th(1,n) = A(pp,1); % Time to reach heliopause [yrs]
            Ve(1,n) = V5(end,1); % Exit velocity [km/s]
            Amax(1,n) = max([max(abs(Ax5)),max(abs(Ay5)),max(abs(Az5))]); % Max acceleration [m/s^2]
            Rmin(1,n) = min(R5); % Minimum solar flyby distance [AU]
        else
            [q, p5] = max(A(:,2));
            te5(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),1)/yr; % elapsed time [yrs] (converted from days)
            R5(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),2)/AU; % orbital distance [AU] (converted from km)
            SMA5(1:min([length(A),length(B),length(C)]),1) = A(1:min([length(A),length(B),length(C)]),3)/AU; % Semi-major axis [AU] (converted from km)
            V5(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),1); % Velocity [km/s]
            Ax5(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),2)*1e3; % Normal acceleration [m/s]
            Ay5(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),3)*1e3; % Lateral acceleration [m/s]
            Az5(1:min([length(A),length(B),length(C)]),1) = B(1:min([length(A),length(B),length(C)]),4)*1e3; % Vertical acceleration [m/s]
            E15(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),1); % Euler angle 1 [deg]
            E25(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),2); % Euler angle 2 [deg]
            E35(1:min([length(A),length(B),length(C)]),1) = C(1:min([length(A),length(B),length(C)]),3); % Euler angle 3 [deg]
            E1r5(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),1); % Euler angle rate 1 [deg]
            E2r5(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),2); % Euler angle rate 2 [deg]
            E3r5(1:min([length(A),length(B),length(C)]),1) = D(1:min([length(A),length(B),length(C)]),3); % Euler angle rate 3 [deg]
            X5(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),1); % X coordinate in Sun-fixed inertial coordinates [km]
            Y5(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),2); % Y coordinate in Sun-fixed inertial coordinates [km]
            Z5(1:min([length(A),length(B),length(C)]),1) = E(1:min([length(A),length(B),length(C)]),3); % Z coordinate in Sun-fixed inertial coordinates [km]
            % Calculate sun-spacecraft angles and instantaneous spacecraft
            % temperatures
            for m = 1:length(X5) % for each data point
                chi25(m,1) = atand(Z5(m,1)/X5(m,1)) + E25(m,1); % inertial plane angle to body lateral plane
                chi35(m,1) = atand(Y5(m,1)/X5(m,1))+E35(m,1); % inertial lateral plane angle to body lateral plane
                % determining whether the sun-facing side is reflective
                % (front) or emmissive (rear)
                if chi25(m,1) < 90 && chi25(m,1) > -90 % front
                    if chi35(m,1) < 90 && chi35(m,1) > -90 % front
                        Cri = 1;
                    else % rear
                        Cri = 2;
                    end
                else % rear
                    if chi35(m,1) < 90 && chi35(m,1) > -90 % front
                        Cri = 2;
                    else % front
                        Cri = 1;
                    end
                end
                % emmissive equilibrium temperature
                T5(m,1) = (((1-Cr(Cri))*I0*abs(sind(chi25(m,1))*sind(chi35(m,1))))/((2-(Cr(1)+Cr(2)))*SBc*R5(m,1)))^0.25;
            end
            
            
            % Overall orbit properties
            Tmax(1,n) = max(T5); % Maximum instantaneous temperature [K]
            th(1,n) = NaN; % Time to reach heliopause [yrs]
            Ve(1,n) = NaN; % Exit velocity [km/s]
            Amax(1,n) = max([max(abs(Ax5)),max(abs(Ay5)),max(abs(Az5))]); % Max acceleration [m/s^2]
            Rmin(1,n) = min(R5) ; % Minimum solar flyby distance [AU]
        end
    end
       
    
    % Update storage counter variable
    n = n + 1;
    %clear a b c d e aID bID cID dID eID A B C D E q m
end

%% Data processing

%% Plots
close all;
%th = A(1,end); % Time to heliopause in days
th = th / yr; % Convert to years
% figure(1)
% hold on;
% box on;
% % grid minor;
% plot(te1(1:500:end,1),R1(1:500:end,1),'r-','LineWidth',2);
% plot(te2(1:500:end,1),R2(1:500:end,1),'k:','LineWidth',2);
% plot(te3(1:500:end,1),R3(1:500:end,1),'b-.','LineWidth',2);
% plot(te4(1:500:end,1),R4(1:500:end,1),'m--','LineWidth',2);
% plot(te5(1:500:end,1),R5(1:500:end,1),'g.','LineWidth',2);
% plot([0,200],[123,123],'k--')
% xlim([0,100]);
% ylim([0,130]);
% set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
% % title('Orbital distance vs. elapsed time','interpreter','latex','fontsize',30);
% ylabel('$R$ [AU]','interpreter','latex','fontsize',45);
% xlabel('$t_e$ [yrs]','interpreter','latex','fontsize',45);
% legend('$m = 2.27$kg','$m = 3.00$kg','$m = 5.00$kg','$m = 10.00$kg','$m = 1.50$kg','Heliopause','interpreter','latex','fontsize',25);
% hold off;
% set(gcf, 'Position', get(0,'Screensize'));
% exportgraphics(gcf,'Mass R vs. t.png','Resolution',200);
% 
% 
% figure(2)
% hold on;
% box on;
% % grid minor;
% plot(te1(1:500:end,1),V1(1:500:end,1),'r-','LineWidth',2);
% plot(te2(1:500:end,1),V2(1:500:end,1),'k:','LineWidth',2);
% plot(te3(1:500:end,1),V3(1:500:end,1),'b-.','LineWidth',2);
% plot(te4(1:500:end,1),V4(1:500:end,1),'m--','LineWidth',2);
% plot(te5(1:500:end,1),V5(1:500:end,1),'g.','LineWidth',2);
% xlim([0,100]);
% set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
% % title('Velocity magnitude vs. elapsed time','interpreter','latex','fontsize',30);
% ylabel('$V$ [kms$^{-1}$]','interpreter','latex','fontsize',45);
% xlabel('$t_e$ [yrs]','interpreter','latex','fontsize',45);
% legend('$m = 2.27$kg','$m = 3.00$kg','$m = 5.00$kg','$m = 10.00$kg','$m = 1.50$kg','interpreter','latex','fontsize',25);
% hold off;
% set(gcf, 'Position', get(0,'Screensize'));
% exportgraphics(gcf,'Mass V vs. t.png','Resolution',200);
% 
% 
% figure(3)
% hold on;
% box on;
% % grid minor;
% plot(te1(1:500:end,1),T1(1:500:end,1),'r-','LineWidth',2);
% plot(te2(1:500:end,1),T2(1:500:end,1),'k:','LineWidth',2);
% plot(te3(1:500:end,1),T3(1:500:end,1),'b-.','LineWidth',2);
% plot(te4(1:500:end,1),T4(1:500:end,1),'m--','LineWidth',2);
% plot(te5(1:500:end,1),T5(1:500:end,1),'g.','LineWidth',2);
% xlim([0,100]);
% set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
% % title('Temperature vs. elapsed time','interpreter','latex','fontsize',30);
% ylabel('$T$ [K]','interpreter','latex','fontsize',45);
% xlabel('$t_e$ [yrs]','interpreter','latex','fontsize',45);
% legend('$m = 2.27$kg','$m = 3.00$kg','$m = 5.00$kg','$m = 10.00$kg','$m = 1.50$kg','interpreter','latex','fontsize',25);
% hold off;
% set(gcf, 'Position', get(0,'Screensize'));
% exportgraphics(gcf,'Mass T vs. t.png','Resolution',200);
% 
% 
% figure(4)
% subplot(2,1,1);
% hold on;
% box on;
% % grid minor;
% plot(te1(1:500:end,1),chi21(1:500:end,1),'r-','LineWidth',2);
% plot(te2(1:500:end,1),chi22(1:500:end,1),'k:','LineWidth',2);
% plot(te3(1:500:end,1),chi23(1:500:end,1),'b-.','LineWidth',2);
% plot(te4(1:500:end,1),chi24(1:500:end,1),'m--','LineWidth',2);
% plot(te5(1:500:end,1),chi25(1:500:end,1),'g.','LineWidth',2);
% xlim([0,100]);
% set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
% % title('$\chi_{2}$ vs. elapsed time','interpreter','latex','fontsize',30);
% ylabel('$\chi_2$ [$^{\circ}$]','interpreter','latex','fontsize',45);
% xlabel('$t_e$ [yrs]','interpreter','latex','fontsize',45);
% legend('$m = 2.27$kg','$m = 3.00$kg','$m = 5.00$kg','$m = 10.00$kg','$m = 1.50$kg','interpreter','latex','fontsize',30);
% hold off;
% 
% subplot(2,1,2);
% hold on;
% box on;
% % grid minor;
% plot(te1(1:500:end,1),chi31(1:500:end,1),'r-','LineWidth',2);
% plot(te2(1:500:end,1),chi32(1:500:end,1),'k:','LineWidth',2);
% plot(te3(1:500:end,1),chi33(1:500:end,1),'b-.','LineWidth',2);
% plot(te4(1:500:end,1),chi34(1:500:end,1),'m--','LineWidth',2);
% plot(te5(1:500:end,1),chi35(1:500:end,1),'g.','LineWidth',2);
% xlim([0,100]);
% set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
% % title('$\chi_{3}$ vs. elapsed time','interpreter','latex','fontsize',30);
% ylabel('$\chi_3$ [$^{\circ}$]','interpreter','latex','fontsize',45);
% xlabel('$t_e$ [yrs]','interpreter','latex','fontsize',45);
% legend('$m = 2.27$kg','$m = 3.00$kg','$m = 5.00$kg','$m = 10.00$kg','$m = 1.50$kg','interpreter','latex','fontsize',30);
% hold off;
% set(gcf, 'Position', get(0,'Screensize'));
% exportgraphics(gcf,'Mass Chi vs. t.png','Resolution',200);
% 
% 
% figure(5)
% subplot(3,1,1);
% hold on;
% box on;
% % grid minor;
% plot(te1(1:500:end,1),Ax1(1:500:end,1),'r-','LineWidth',2);
% plot(te2(1:500:end,1),Ax2(1:500:end,1),'k:','LineWidth',2);
% plot(te3(1:500:end,1),Ax3(1:500:end,1),'b-.','LineWidth',2);
% plot(te4(1:500:end,1),Ax4(1:500:end,1),'m--','LineWidth',2);
% plot(te5(1:500:end,1),Ax5(1:500:end,1),'g.','LineWidth',2);
% xlim([0,100]);
% set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
% % title('Acceleration in the body-fixed normal axis x vs. elapsed time','interpreter','latex','fontsize',30);
% ylabel('$A_x$ [ms$^{-2}$]','interpreter','latex','fontsize',45);
% xlabel('$t_e$ [yrs]','interpreter','latex','fontsize',45);
% legend('$m = 2.27$kg','$m = 3.00$kg','$m = 5.00$kg','$m = 10.00$kg','$m = 1.50$kg','interpreter','latex','fontsize',30);
% hold off;
% 
% subplot(3,1,2);
% hold on;
% box on;
% % grid minor;
% plot(te1(1:500:end,1),Ay1(1:500:end,1),'r-','LineWidth',2);
% plot(te2(1:500:end,1),Ay2(1:500:end,1),'k:','LineWidth',2);
% plot(te3(1:500:end,1),Ay3(1:500:end,1),'b-.','LineWidth',2);
% plot(te4(1:500:end,1),Ay4(1:500:end,1),'m--','LineWidth',2);
% plot(te5(1:500:end,1),Ay5(1:500:end,1),'g.','LineWidth',2);
% xlim([0,100]);
% set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
% % title('Acceleration in the body-fixed lateral axis y vs. elapsed time','interpreter','latex','fontsize',30);
% ylabel('$A_y$ [ms$^{-2}$]','interpreter','latex','fontsize',45);
% xlabel('$t_e$ [yrs]','interpreter','latex','fontsize',45);
% legend('$m = 2.27$kg','$m = 3.00$kg','$m = 5.00$kg','$m = 10.00$kg','$m = 1.50$kg','interpreter','latex','fontsize',30);
% hold off;
% 
% subplot(3,1,3);
% hold on;
% box on;
% % grid minor;
% plot(te1(1:500:end,1),Az1(1:500:end,1),'r-','LineWidth',2);
% plot(te2(1:500:end,1),Az2(1:500:end,1),'k:','LineWidth',2);
% plot(te3(1:500:end,1),Az3(1:500:end,1),'b-.','LineWidth',2);
% plot(te4(1:500:end,1),Az4(1:500:end,1),'m--','LineWidth',2);
% plot(te5(1:500:end,1),Az5(1:500:end,1),'g.','LineWidth',2);
% xlim([0,100]);
% set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
% % title('Acceleration in the body-fixed vertical axis z vs. elapsed time','interpreter','latex','fontsize',30);
% ylabel('$A_z$ [ms$^{-2}$]','interpreter','latex','fontsize',45);
% xlabel('$t_e$ [yrs]','interpreter','latex','fontsize',45);
% legend('$m = 2.27$kg','$m = 3.00$kg','$m = 5.00$kg','$m = 10.00$kg','$m = 1.50$kg','interpreter','latex','fontsize',30);
% hold off;
% set(gcf, 'Position', get(0,'Screensize'));
% exportgraphics(gcf,'Mass A vs. t.png','Resolution',200);
% 
% 
% figure(6)
% subplot(1,2,1)
% hold on;
% box on;
% % grid minor;
% plot(Y1(1:500:end,1)/AU,Z1(1:500:end,1)/AU,'r-','LineWidth',2);
% plot(Y2(1:500:end,1)/AU,Z2(1:500:end,1)/AU,'k:','LineWidth',2);
% plot(Y3(1:500:end,1)/AU,Z3(1:500:end,1)/AU,'b-.','LineWidth',2);
% plot(Y4(1:500:end,1)/AU,Z4(1:500:end,1)/AU,'m--','LineWidth',2);
% plot(Y5(1:500:end,1)/AU,Z5(1:500:end,1)/AU,'g.','LineWidth',2);
% set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
% % title('ZY-plane trajectory','interpreter','latex','fontsize',30);
% ylabel('$Z$ [AU]','interpreter','latex','fontsize',45);
% xlabel('$Y$ [AU]','interpreter','latex','fontsize',45);
% legend('$m = 2.27$kg','$m = 3.00$kg','$m = 5.00$kg','$m = 10.00$kg','$m = 1.50$kg','interpreter','latex','fontsize',30);
% hold off;
% 
% 
% 
% subplot(1,2,2)
% hold on;
% box on;
% % grid minor;
% plot(X1(1:500:end,1)/AU,Y1(1:500:end,1)/AU,'r-','LineWidth',2);
% plot(X2(1:500:end,1)/AU,Y2(1:500:end,1)/AU,'k:','LineWidth',2);
% plot(X3(1:500:end,1)/AU,Y3(1:500:end,1)/AU,'b-.','LineWidth',2);
% plot(X4(1:500:end,1)/AU,Y4(1:500:end,1)/AU,'m--','LineWidth',2);
% plot(X5(1:500:end,1)/AU,Y5(1:500:end,1)/AU,'g.','LineWidth',2);
% set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
% % title('XY-plane trajectory','interpreter','latex','fontsize',30);
% ylabel('$Y$ [AU]','interpreter','latex','fontsize',45);
% xlabel('$X$ [AU]','interpreter','latex','fontsize',45);
% legend('$m = 2.27$kg','$m = 3.00$kg','$m = 5.00$kg','$m = 10.00$kg','$m = 1.50$kg','interpreter','latex','fontsize',30);
% hold off;
% set(gcf, 'Position', get(0,'Screensize'));
% exportgraphics(gcf,'Mass 3D trajectory.png','Resolution',200);
% 

figure(1)

subplot(3,2,1);
hold on;
box on;
% grid minor;
plot(mass(1:CurrentSimNo),th,'rv','LineWidth',3,'MarkerSize',15);
xlim([min(mass)/1.2,max(mass)*1.03]);
ylim([min(th)/2,max(th)*1.075]);
set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
% title('Solar exit time vs. spacecraft mass','interpreter','latex','fontsize',30);
ylabel('$t_h$ [yrs]','interpreter','latex','fontsize',40);
xlabel('$m_{sc}$ [kg]','interpreter','latex','fontsize',40);
%legend('Measured data','Second-order fit','interpreter','latex','fontsize',30);
hold off;


subplot(3,2,2);
hold on;
box on;
% grid minor;
plot(mass(1:CurrentSimNo),Rmin,'b^','LineWidth',3,'MarkerSize',15);
set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
xlim([min(mass)/1.13,max(mass)*1.03]);
ylim([min(Rmin)/1.16,max(Rmin)*1.075]);
% title('Minimum solar flyby distance vs. spacecraft mass','interpreter','latex','fontsize',30);
ylabel('$R_{min}$ [AU]','interpreter','latex','fontsize',40);
xlabel('$m_{sc}$ [kg]','interpreter','latex','fontsize',40);
%legend('Measured data','Second-order fit','interpreter','latex','fontsize',30);
hold off;

subplot(3,2,3);
hold on;
box on;
% grid minor;
plot(mass,Tmax,'k>','LineWidth',3,'MarkerSize',15,'HandleVisibility','off');
plot([min(mass),max(mass)],[443,443],'k--','LineWidth',3);
set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
xlim([min(mass)/1.13,max(mass)*1.03]);
ylim([min(Tmax)/1.16,max(Tmax)*1.075]);
% title('Maximum temperature vs. spacecraft mass','interpreter','latex','fontsize',30);
ylabel('$T_{s_{max}}$ [K]','interpreter','latex','fontsize',40);
xlabel('$m_{sc}$ [kg]','interpreter','latex','fontsize',40);
lt = legend('$T_{op_{max}}$','interpreter','latex','fontsize',30);
lt.LineWidth = 2;
hold off;

subplot(3,2,4);
hold on;
box on;
% grid minor;

plot(mass,Amax*1e3,'m<','LineWidth',3,'MarkerSize',15);
xlim([min(mass)/1.13,max(mass)*1.03]);
ylim([min(Amax*1e3)/1.16,max(Amax*1e3)*1.075]);
set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
% title('Maximum acceleration vs. spacecraft mass','interpreter','latex','fontsize',30);
ylabel('$A_{max}$ [mms$^{-2}$]','interpreter','latex','fontsize',30);
xlabel('$m_{sc}$ [kg]','interpreter','latex','fontsize',30);
%legend('Measured data','Second-order fit','interpreter','latex','fontsize',30);
hold off;

subplot(3,2,5);
hold on;
box on;
% grid minor;
plot(mass,Ve,'go','LineWidth',3,'MarkerSize',15);
set(gca,'LineWidth',4.5,'FontSize',30,'TickLabelInterpreter','LaTeX');
xlim([min(mass)/1.13,max(mass)*1.03]);
ylim([min(Ve)/1.16,max(Ve)*1.075]);
% title('Exit velocity vs. spacecraft mass','interpreter','latex','fontsize',30);
ylabel('$V_{e}$ [kms$^{-1}$]','interpreter','latex','fontsize',40);
xlabel('$m_{sc}$ [kg]','interpreter','latex','fontsize',40);
%legend('Measured data','Second-order fit','interpreter','latex','fontsize',30);
hold off;
set(gcf, 'Position', get(0,'Screensize'));
exportgraphics(gcf,'Mass sensitivity.png','Resolution',200);


% figure(8)
% hold on;
% box on;
% % grid minor;
% plot(te1(1:500:end,1),log10(SMA1(1:500:end,1)),'r-','LineWidth',2);
% plot(te2(1:500:end,1),log10(SMA2(1:500:end,1)),'k:','LineWidth',2);
% plot(te3(1:500:end,1),log10(SMA3(1:500:end,1)),'b-.','LineWidth',2);
% plot(te4(1:500:end,1),log10(SMA4(1:500:end,1)),'m--','LineWidth',2);
% plot(te5(1:500:end,1),log10(SMA5(1:500:end,1)),'g.','LineWidth',2);
% xlim([0,100]);
% % title('Semi-major axis vs. elapsed time','interpreter','latex','fontsize',30);
% ylabel('log$_{10}$(SMA) [AU]','interpreter','latex','fontsize',30);
% xlabel('$t_e$ [yrs]','interpreter','latex','fontsize',30);
% legend('$m = 2.27$kg','$m = 3.00$kg','$m = 5.00$kg','$m = 10.00$kg','$m = 1.50$kg','interpreter','latex','fontsize',30);
% hold off;
% set(gcf, 'Position', get(0,'Screensize'));
% exportgraphics(gcf,'Mass log10(SMA).png','Resolution',200);


% %% Output file
% 
% ODID = fopen('OrbitalDataV5.txt','w');
% fprintf(ODID,'%19s %24s %37s %24s\n','Spacecraft mass [kg]','Time to Heliopause [yrs]','Minimum solar approach [AU]','Maximum temperature [K]');
% 
% for j = 1:CurrentSimNo
%     if ~isnan(th(1,j)) && th(1,j) < 100
%         fprintf(ODID,'%17.2f %24d %34.3f %24d\n',mass(j),floor(th(1,j)),Rmin(j),floor(Tmax(j)));
%     else
%         fprintf(ODID,'%17.2f %24s %34.3f %24d\n',mass(j),' Not extrasolar',Rmin(j),floor(Tmax(j)));
%     end
% end
% fclose(ODID);





