%% Define constants
Lv = 2.5e6; % [J / kg]
densityWater = 1e3; % [kg / m^3]
% densityAir = 1.2922; % [kg / m^3] (Wikipedia)
% densityAir = 1.2; % [kg / m^3] (http://www.engineeringtoolbox.com/air-density-specific-weight-d_600.html, approx. at 20 C [68 F])
densityAir = 1.28; % WRF
% cpAir = 1005; % [J / kg K] (http://www.engineeringtoolbox.com/air-properties-d_156.html at 20 C)
cpAir = 1004.6; % WRF
% R = 287; % [J / kg K]
R = 287.04; % [J / kg K] WRF
% Rv = 461; % [J / kg K]
Rv = 461.6; % [J / kg K] WRF
g = 9.81; % [m / s^2]

% Net LW calculation constants (after Stevens, et al., 2005 - Evaluation of large-Eddy simulations...)
kappa = 85; % [m^2 / kg]
alpha = 1; % [m^(-4/3)]
D = 3.75 * 10^-6; % [1 / s]

plotDir = '/home/hyang/Dropbox/thesisPlots_v2/';
png.format = 'png';

plotList = {'r', 'b--', 'g', 'o', 'm'};