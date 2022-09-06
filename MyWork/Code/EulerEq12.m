% Function that provides angular acceleration as a result of a specific
% torque.
% Tanmay Ubgade | 220531

function [omega_dot] = EulerEq12(T,I_tot)
%% Inputs
% Torque vector where T(1) is T1 (x-ax torque) and T(2) is T2 (y-ax torque)
% I is the 2nd moment of inertia of the sailcraft, acquired from the
% MomofInertia function

%% Processing
I = diag(I_tot);
It = I(1);
omega_dot = T/It;
end