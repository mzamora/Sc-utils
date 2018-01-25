i = 1; j = 1; t = 1;

%% input_sounding
fid = fopen('sound_in', 'w');

fprintf( fid, '%.2f, %.2f, %.2f, %.2f, %.2f\n', wrf.PSFC( i, j, t ) ./ 100, wrf.T2( i, j, t ), wrf.Q2( i, j, t ).*1000, wrf.U10( i, j, t ), wrf.V10( i, j, t ) );

for z = 1:size(wrf.T, 3)
	fprintf( fid, '%.2f, %.2f, %.2f, %.2f, %.2f\n', wrf.height( i, j, z, t ), wrf.T( i, j, z, t )+300, (wrf.QVAPOR( i, j, z, t ) + wrf.QCLOUD( i, j, z, t ) + wrf.QRAIN( i, j, z, t) ).*1000, wrf.U( i, j, z, t ), wrf.V( i, j, z, t ) );
end
fclose(fid);
