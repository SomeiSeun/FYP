% Calcuating the moment of interia of the Solar Sailcraft
% By Tanmay Ubgade 220529

function I_tot = MomOfInertia()

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

    I_tot = (I_Sail + I_CubeSat)*(1000^2)/1000; % In kg/m^2
end