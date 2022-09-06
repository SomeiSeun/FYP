% Script to calculate the total power requirements for sailcraft
% Tanmay Ubgade 220321

%% housekeeping
clear
clc
close all

%% Constants

k = 1.38*10^-20; % Boltzmann constant
T_S = 28.5; % System equivalent noise temperature in K
N0 = k*T_S; % Noise spectral density
c = 299792458; % Speed of light
AU = 1495978.70691; % AU to km conversion

%% Base parameters

% Signal
%lambda = 3.75e-2; % Wavelength in m (X-band is 3.75 - 2.4 cm)
freq = 8.4*10^9; % 8.4 GHZ 
lambda = c/freq;
% Distance
r = 7*AU*10^3; % distance of sailcraft from Earth
A = 4*pi*r^2; % Area of sphere

% Transmitter

P_T = 4; % Power transmitted
%G_T = 4*pi*A_T/lambda^2; % Gain of transmitting antenna
G_T = 6.5*10^4; % Voyager gain at Jupiter
A_T = G_T*(lambda^2)/(4*pi);

% Receiver

rho = P_T*G_T/A; % Power flux density at receiver
mu = 0.9; % Receiving antenna efficiency
%A_R = 1; % Area of receiving antenna 
%P_R = rho*mu*A_R; % Received signal
%G_R = 4*pi*A_R/lambda^2; % Receiving antenna gain
SNR = 5.48*10^5; %SNR for voyager at 8.4GHZ, SNR = P_R/(L * N0)
L = 0.7; % Approximate total loss
P_R = SNR*L*N0;
A_R = P_R/rho;
G_R = 4*pi*A_R/lambda^2;

RT_ratio = G_T*G_R*(lambda^2)/(4*pi*r)^2;
P_TTC = P_R/RT_ratio;

%% C&DH

P_CDH = 0.4;

%% Link Design Equation
% EIRP = P + L_l + G_t

P_db = 10*log10(4);
L_l = -1; % in dB
G_t = 44.3 - 10*log10(74^2);
EIRP = P_db + L_l + G_t;

% Changing data rate
%datarate = 1:1:2048; % in bps

% Changing range 
r = linspace(0,8)*AU*1000;
datarate = 1024;


L_s = 147.55 - 20*log10(r)-20*log10(freq);
L_alp = 4*10^-2; % in dB taken from fig 13.10 in SMAD
G_r = -159.59 + 20*log10(35) + 20*log10(freq) + 10*log(0.663); % DSN antenna radius is 35m
T_s = 135; % 135 downlink, 614 uplink



EbNo = EIRP + L_s + L_alp + G_r + 228.6 - 10*log10(T_s) - 10*log10(datarate);

%{
fig1 = figure(1);
hold on 
plot(datarate,EbNo)
plot([1024,1024],[0,45],'--')
xlabel('Data Rate [bps]')
ylabel('{E_b}/{N_o}')
box on
grid on
grid minor
ylim([min(EbNo) 35])
xticks([0 512 1024 1024+512 2048])
hold off

fig1.Units = 'inches';
fig1.Position(3) = 2.8;
fig1.Position(4) = 2.8;
set(fig1.Children, 'FontName', 'Arial', 'FontSize', 11);
legend('{E_b}/{N_o}','1024 bps','Location','Northwest','FontSize', 8)
%print('LinkDesignDataRate', '-depsc')
%}

%
fig1 = figure(1);
hold on
plot(r/(1000*AU),EbNo)
plot([7,7],[0,55],'--')
hold off
xlabel('Range [AU]')
ylabel('{E_b}/{N_o}')
xlim([0 8])
ylim([0 50])
box on
grid on
grid minor
fig1.Units = 'inches';
fig1.Position(3) = 2.8;
fig1.Position(4) = 2.8;
set(fig1.Children, 'FontName', 'Arial', 'FontSize', 11);
legend('{E_b}/{N_o}','7 AU limit','Location','Northwest','FontSize', 8)
print('LinkDesignRange', '-depsc')

%}


% boltzmann constant in dB scale = 228.6

%% Solar Panel Sizing

% CubeSat Panel assumptions
P_req = 5; % in W
P_sa = P_req/0.8;

% 19-20 years in Space

% Efficiencies: Si = 14.8%, GaAs = 18.5%, Multijunction = 22%
I_Earth = 1367; %W/m^2
I_limit = I_Earth/(2^2); % intensity at 8 AU, inverse square law

P0_SI = 0.148*I_limit;
P0_GaAs = 0.185*I_limit;
P0_mult = 0.22*I_limit;
P0 = [P0_SI, P0_GaAs,P0_mult]
P_BOL = P0*0.77*cosd(60); % Inherent degredation = 0.77
P_EOL = P_BOL.*(1 - [0.0375, 0.0275, 0.005]).^20;

area_sa = 1./(P_EOL./P_sa)
mass_sa1 = 0.04*P0