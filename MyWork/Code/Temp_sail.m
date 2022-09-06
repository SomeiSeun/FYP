% Function that calculates temperature of the Sail based on emmisivity and
% orientation.

function T = Temp_sail(r,chi2,chi3,rho_r,rho_e)
%% Inputs
% r - Radius vector of Sun-Sail positions
% Chi2 -
% Chi3 -
% rho_r - reflectivity of reflective surface
% rho_r - reflectivity of emissive surface
SBc = 5.670374419e-8; % Stefan-Boltzmann constant

%% Calculations

T = ((1-rho_r)*SolarIntensity(r).*abs(sind(chi2).*sind(chi3))/((2-rho_r-rho_e)*SBc)).^(1/4);
end