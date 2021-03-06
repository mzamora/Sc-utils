%% Computes potential temperature given temperature, pressure, and reference pressure (100 kPa)

function out = theta( T, P )

%% Define constants
cpAir = 1005; % [J / kg K] (http://www.engineeringtoolbox.com/air-properties-d_156.html at 20 C)
R = 287; % [J / kg K]
P0 = 100000; % [Pa]
%% Go

out = T .* ( P0 ./ P ) .^ (R ./ cpAir);

end