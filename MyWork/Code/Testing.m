% housekeeping
clear all
close all
clc

    i = 1;    
    n = 1;
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
    fscanf(aID,'%s %s %s',[3,1]);
    fscanf(bID,'%s %s %s %s',[4,1]);
    fscanf(cID,'%s %s %s',[3,1]);
    fscanf(dID,'%s %s %s',[3,1]);
    fscanf(eID,'%s %s %s',[3,1]);
    
    % Store data into file-specific matrices
    A = transpose(fscanf(aID,'%f %f %f',[3,inf]));
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

    te = A(1:end,1); % elapsed time [yrs] (converted from days)
    R = A(1:end,2); % orbital distance [AU] (converted from km)
    SMA = A(1:end,3); % Semi-major axis [AU] (converted from km)
    V = B(1:end,1); % Velocity [km/s]
    Ax = B(1:end,2)*1e3; % Normal acceleration [m/s]
    Ay = B(1:end,3)*1e3; % Lateral acceleration [m/s]
    Az = B(1:end,4)*1e3; % Vertical acceleration [m/s]
    E1 = C(1:end,1); % Euler angle 1 [deg]
    E2 = C(1:end,2); % Euler angle 2 [deg]
    E3 = C(1:end,3); % Euler angle 3 [deg]
    E1_dot = D(1:end,1); % Euler angle rate 1 [deg]
    E2_dot = D(1:end,2); % Euler angle rate 2 [deg]
    E3_dot = D(1:end,3); % Euler angle rate 3 [deg]
    X = E(1:end,1); % X coordinate in Sun-fixed inertial coordinates [km]
    Y = E(1:end,2); % Y coordinate in Sun-fixed inertial coordinates [km]
    Z = E(1:end,3); % Z coordinate in Sun-fixed inertial coordinates [km]
    % Calculate sun-spacecraft angles and instantaneous spacecraft
    % temperatures
    for m = 1:length(X) % for each data point
        chi2(m,1) = atand(Z(m,1)/X(m,1)) + E2(m,1); % inertial plane angle to body lateral plane
        chi3(m,1) = atand(Y(m,1)/X(m,1))+E3(m,1); % inertial lateral plane angle to body lateral plane
        % determining whether the sun-facing side is reflective
        % (front) or emmissive (rear)
        if chi2(m,1) < 90 && chi2(m,1) > -90 % front
            if chi3(m,1) < 90 && chi3(m,1) > -90 % front
                Cri = 1;
            else % rear
                Cri = 2;
            end
        else % rear
            if chi3(m,1) < 90 && chi3(m,1) > -90 % front
                Cri = 2;
            else % front
                Cri = 1;
            end
        end
        % emmissive equilibrium temperature
        Temp(m,1) = (((1-Cr(Cri))*I0*abs(sind(chi2(m,1))*sind(chi3(m,1))))/((2-(Cr(1)+Cr(2)))*SBc*R(m,1)))^0.25;
    end

    % Overall orbit properties
    Tmax(1,n) = max(Temp); % Maximum instantaneous temperature [K]
    [qq, pp] = max(A(:,2));
    th = A(pp,1); % Time to reach heliopause [yrs]
    Ve = V(end,1); % Exit velocity [km/s]
    Amax = max([max(abs(Ax)),max(abs(Ay)),max(abs(Az))]); % Max acceleration [m/s^2]
    Rmin = min(R); % Minimum solar flyby distance [AU]
    
    figure(1)
    plot(te,R)
    xlabel('Time in seconds')
    ylabel('Radius in km')