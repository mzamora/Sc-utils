%% Generate WRF input_sounding from loaded UCLALES data
% Assume data has already been loaded (use go.m)

%% Constants
Lv = 2.5e6; % [J / kg]
cpAir = 1005; % [J / kg K] (http://www.engineeringtoolbox.com/air-properties-d_156.html at 20 C)
R = 287; % [J / kg K]

%% Setup

% Define timestep to retrieve data from
t = 1;

% type = 'qv';
type = 'qt';

% Define output directory
outputDir = ['scm/uclales/' type '/']; if ~exist(outputDir, 'dir'), mkdir(outputDir); end

%% Compute potential temperature theta from liquid water potential temperature theta_l

for z = 1:size(ps.t, 1)
	ps.theta(z, t) = ps.t(z, t) + (Lv ./ cpAir) .* ( ( ps.p(1, t) ./ ps.p(z, t) ) .^ (R ./ cpAir) ) .* (ps.l(z, t) ./ 1000);
end

%% Make sounding

%% Make input_sounding for WRF SCM

% input_sounding
fid = fopen([outputDir '/input_sounding'], 'w');
% z_surface, u10, v10, temp2, q2, psfc
fprintf( fid, '0.0, %.1f, %.1f, %.1f, %.6f, %.1f\n', ps.u(1, t), ps.v( 1, t ), ps.theta( 1, t ), ps.q( 1, t ) ./ 1000, ps.p( 1, t ));

for z = 2:size(ps.zt, 1)
	% height, u, v, theta, qv
	if strcmpi( type, 'qt' )
		fprintf( fid, '%.1f, %.1f, %.1f, %.1f, %.6f\n', ps.zt( z, t ), ps.u( z, t ), ps.v( z, t ), ps.theta( z, t ), ps.q( z, t ) ./ 1000 );
	elseif strcmpi( type, 'qv' )
		fprintf( fid, '%.1f, %.1f, %.1f, %.1f, %.6f\n', ps.zt( z, t ), ps.u( z, t ), ps.v( z, t ), ps.theta( z, t ), ( ps.q( z, t ) - ps.l(z, t) ) ./ 1000 );
	else
		error('Undefined sounding type (qt or qv)')
	end
end
fclose(fid);