%% Reads in UCLALES sounding (sound_in) and converts to WRF-SCM input sounding (input_sounding)

% Path to UCLALES sounding
input_sounding = '/home/hyang/research/matlab/data/LES/sound_in_minus1';

% Path to output sounding for WRF-SCM
output_sounding = '/home/hyang/research/matlab/projects/scm_qt_qv_comparison/sounding_CGILS_qt_minus1';

%% Read in UCLALES sound_in
fid = fopen(input_sounding, 'r');

% P_sfc			Theta_l		qt(g/kg)	u	v
% height(m)		Theta_l		qt(g/kg)	u	v
uclales = textscan(fid, '%f %f %f %f %f');

height = uclales{1}; theta_l = uclales{2}; qt = uclales{3} ./ 1000; u = uclales{4}; v = uclales{5};

% Theta_l = theta in the absence of ql

fclose(fid);

%% Make WRF-SCM input_sounding

% input_sounding
fid = fopen(output_sounding, 'w');
fprintf( fid, '0.0, %.1f, %.1f, %.1f, %.6f, %.1f\n', u(1), v(1), theta_l(1), qt(1), height(1).*100);
for z = 2:length(height)
	fprintf( fid, '%.1f, %.1f, %.1f, %.1f, %.6f\n', height(z), u(z), v(z), theta_l(z), qt(z));
end
fclose(fid);