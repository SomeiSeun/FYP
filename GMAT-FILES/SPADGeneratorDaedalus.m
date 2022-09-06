% This code estimates solar radiation pressure parameters for the Sunbeam
% spacecraft and tabulates it with azimuth and elevation (orbital distance)
% wrt the Sun
% Written by Filippos Geragidis on 14/2/2021

% Housekeeping
clear;
clc;

% Definition of sailcraft properties
Cr = [0.91,0.16]; % reflective and emmissive side coefficients of reflectivity
ASRP = 454; % solar radiation pressure area in m^2

% Initialisations
n = 1; % data storage index

% Data generation loop
for a = -180:5:180 % for azimuth from -180deg to +180deg
    
    for e = -90:5:90 % for elevation from -90deg to +90deg
        
        % reflectivity assignment
        if e >=0 % front side
            Cri = 1;
        else % back side
            Cri = 2;
        end
        
        % Estimation of total force scale factor (m^2)
        Fsolar = (1+Cr(Cri)) * ASRP * sind(e);
        % Transformation of total force vector into components and storage
        Fy(n,1) = 0;
        Fz(n,1) = 0;
        Fx(n,1) = -Fsolar * sind(e);
        
        % Storage of azimuth and elevation
        A(n,1) = a;
        E(n,1) = e;
        
        % Index tick-over
        n = n + 1;
        
    end
    
end

% Creation of total matrix for inspection
SPAD_matrix= [A, E, Fx, Fy, Fz];

% SunbeamSPAD.dat file creation
filename='SunbeamSPADV8c.dat';
SPADID = fopen(filename,'w');
fprintf(SPADID,'%-20s %2s %-20s \n','Version',':','1.0');
fprintf(SPADID,'%-20s %2s %-20s \n','System',':','Sunbeam');
fprintf(SPADID,'%-20s %2s %-20s \n','Analysis Type',':','Area');
fprintf(SPADID,'%-20s %2s %-20s \n','Current time',':',datetime);
fprintf(SPADID,'\n');
fprintf(SPADID,'\n');
fprintf(SPADID,'%-10s %2s %-1d \n','Motion',':',1);
fprintf(SPADID,'  %-8s %2s %-7s \n','Name',':','Azimuth');
fprintf(SPADID,'  %-8s %2s %-7s \n','Method',':','Step');
fprintf(SPADID,'  %-8s %2s %-3d \n','Minimum',':',-180);
fprintf(SPADID,'  %-8s %2s %-3d \n','Maximum',':',180);
fprintf(SPADID,'  %-8s %2s %-1d \n','Step',':',5);
fprintf(SPADID,'%-10s %2s %-1d \n','Motion',':',2);
fprintf(SPADID,'  %-8s %2s %-7s \n','Name',':','Elevation');
fprintf(SPADID,'  %-8s %2s %-7s \n','Method',':','Step');
fprintf(SPADID,'  %-8s %2s %-3d \n','Minimum',':',-90);
fprintf(SPADID,'  %-8s %2s %-3d \n','Maximum',':',90);
fprintf(SPADID,'  %-8s %2s %-1d \n','Step',':',5);
fprintf(SPADID,'%5s \n',': END');
fprintf(SPADID,'\n');
fprintf(SPADID,'%12s %1s %5d \n','Record count',':',length(SPAD_matrix));
fprintf(SPADID,'%-10s %-7s %-18s %-18s %-18s \n','Azimuth','Elevation','Force(X)','Force(Y)','Force(Z)');
fprintf(SPADID,'%10s %7s %18s %18s %18s \n','degrees','degrees','m^2','m^2','m^2');
fprintf(SPADID,'%-10s %-7s %-18s %-18s %-18s \n','----------','-------','------------------','------------------','------------------');
fprintf(SPADID,'%10.2f %7.2f %18.14f %18.14f %18.14f \n',SPAD_matrix.');
fclose(SPADID);





