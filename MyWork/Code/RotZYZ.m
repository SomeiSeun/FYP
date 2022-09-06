% Function to provide output of overall rotation matrix (Direction Cosine
% Matrix)
% Tanmay Ubgade | 220531

function [R] = RotZYZ(phi,gam,psi)
%% Inputs
% phi - First rotation about Z axis
% gam - Second rotation about Y axis
% psi - Third rotation about Z axis

%% Rotation matrices
R3_phi = [cos(phi) sin(phi) 0;
         -sin(phi) cos(phi) 0;
          0        0        1];
      
R2_gam = [cos(gam) 0        sin(gam);
          0        1        0;
         -sin(gam) 0        cos(gam)];

R3_psi = [cos(psi) sin(psi) 0;
         -sin(psi) cos(psi) 0;
          0        0        1];

R = R3_phi*R2_gam*R3_psi;

end