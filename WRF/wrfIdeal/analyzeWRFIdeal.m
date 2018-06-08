%% Compare WRF ideal simulations

function cv = analyzeWRFIdeal(filePath, runName)

% stupid hack to make this function work
path{1} = filePath;
name{1} = runName;

% clear all
clc
close all

%% Define physics constants
defineConstants;

%% Define plotDir

plotDir = ['/home/hyang/Dropbox/matlab/projects/wrfIdeal/plots_new/' name{1} '/'];
if ~exist(plotDir), mkdir(plotDir), end

png.format = 'png';
eps.format = 'eps';

%% BLOCK

vars = {'Times', 'PBLH', 'QVAPOR', 'QCLOUD', 'QRAIN', 'RTHRATEN', 'T', 'PH', 'PHB', 'EXCH_H', 'U', 'V', 'W', 'P', 'PB', 'DNW', 'LH', 'HFX', 'HGT', 'ZNW', 'PSFC', 'RAINC', 'RAINNC', 'TSK', ...
	'Q2','T2','TH2', 'MU', 'MUB', 'ZNU', 'ZNW', 'RH'}; 


varsMYNN = {'QKE', 'EL_PBL', 'QSHEAR', 'QBUOY', 'QDISS', 'QWT', 'SH3D', 'DTKE'};
% varsMYNN = [];

for n = 1%:length(path)
	
	for v = 1:length(vars)
		[wrf(n).(vars{v})] = wrfLoad( path{n}, vars(v) );
	end
	if strcmpi(name{n}(1:4), 'mynn')
		for v = 1:length(varsMYNN)
			[wrf(n).(varsMYNN{v})] = wrfLoad( path{n}, varsMYNN(v) );
		end
	end
% 	[wrf(n).Times, wrf(n).PBLH, wrf(n).QVAPOR, wrf(n).QCLOUD, wrf(n).T, wrf(n).PH, wrf(n).PHB, wrf(n).EXCH_H, wrf(n).U, wrf(n).V, ...
% 		wrf(n).P, wrf(n).PB] = wrfLoad( path{n}, vars );
% 	if strcmpi(name{n}(1:4), 'mynn')
% 		[wrf(n).QKE, wrf(n).EL_PBL, wrf(n).QSHEAR, wrf(n).QBUOY, wrf(n).QDISS, wrf(n).QWT, wrf(n).SH3D, wrf(n).DTKE] = wrfLoad( path{n}, varsMYNN );
% 	end
end

%% Fix for registry thing
if ~isfield(wrf, 'T_PHY')
	for i = 1:size(wrf.T, 1)
		for j = 1:size(wrf.T, 2)
			for k = 1:size(wrf.T, 3)
				for t = 1:size(wrf.T, 4)
					wrf.T_PHY(i, j, k, t) = tFromTheta( wrf.T(i, j, k, t) + 300, wrf.P(i, j, k, t) + wrf.PB(i, j, k, t) );
				end
			end
		end
	end
end


disp('Success')
LESps = '~/Dropbox/UCSD/WRF PBL Parameterization/Code & Data/Data/DYCOMS_ocean_run/cgils_s12_ctl.ps.nc';
LESts = '~/Dropbox/UCSD/WRF PBL Parameterization/Code & Data/Data/DYCOMS_ocean_run/cgils_s12_ctl.ts.nc';

%% Load LES
loadLES;

%% Get height & times
% wrf(1).PBLH = wrf(2).PBLH; % have to put this stupid hack or normalization doesn't work for wrong pblh case
for n = 1:length(wrf)
% 	wrf(n).T = wrf(n).T + 300;
	wrf(n).time = datenum(wrf(n).Times');
	
	temp=(wrf(n).PH+wrf(n).PHB)./9.81;
	wrf(n).height=(temp(:,:,1:size(temp,3)-1,:)+temp(:,:,2:size(temp,3),:))/2;
	wrf(n).heightstag = temp;
	
	% Normalize height
% 	for i = 1:size(wrf(1).T, 1)
% 		for j = 1:size(wrf(1).T, 2)
% 			
% % 			for t = 1:size(wrf(1).T, 4)
% % 				wrf(n).heightstag(i,j,1,t) = 0.5 * ( wrf(n).height(i,j,1,t) );
% % 				for k = 2:size(wrf(1).T, 3)
% % 					wrf(n).heightstag(i,j,k,t) = 0.5 * ( wrf(n).height(i,j,k,t) + wrf(n).height(i,j,k-1,t) );
% % 				end
% % 				wrf(n).heightstag(i,j,size(wrf(1).T, 3) + 1,t) = 99999;
% % 			end
% 			
% % 			for k = 1:size(wrf(n).T, 3)
% % 				for t = 2:size(wrf(n).T, 4)
% % 					wrf(n).znorm(i,j,k,t) = wrf(n).height(i,j,k,t) ./ wrf(n).PBLH(i, j, t);
% % 					wrf(n).znormstag(i,j,k,t) = wrf(n).heightstag(i,j,k,t) ./ wrf q(n).PBLH(i, j, t);
% % 				end
% % 			end
% 		end
% 	end
end

%% MLM CV
x = 10:41; y = 10:41;
[cv] = wrfMLCVavg(wrf(1), wrf(1).time(1), wrf(1).time(end), x, y);

%% CV Plot 1
wrf(1).XLONG = wrfLoad(path{1}, {'XLONG'});
wrf(1).XLAT = wrfLoad(path{1}, {'XLAT'});

sqPlotHeight = 0.3;
longPlotWidth = 0.28;
offset = 0.05;

timeMarks = ((1:16:length(cv.time)) - 1) ./ 4;
timeVec = (1:length(cv.time))./4;

figure('position', [120 20 1500 900])

for t = 1:8:length(cv.time)
	clf
	% Left top
	axes('position', [0.04 0.655 0.2 sqPlotHeight])
	plot(timeVec,cv.avg.lwp)
	hold on; plot(les.timeTS./3600, les.lwp./1000, 'k')
	
	xlim([timeMarks(1), timeMarks(end)])
	set(gca,'xtick',timeMarks, 'xticklabel', [])
	grid on
	title('<LWP [kg/m^2]>')
	
	% Left mid
	axes('position', [0.04 0.31 0.2 sqPlotHeight])
	plot(timeVec,cv.avg.PBLH)
	hold on; plot(les.timeTS./3600, les.zi1, 'k')
	
	xlim([timeMarks(1), timeMarks(end)])
	set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
	grid on
	title('<PBLH [m]>')
	
	% Center LWP
	axes('position', [0.28 0.31 0.45 0.645])
	n = 1;
	for i = 1:size(wrf(1).QCLOUD, 1) % Compute LWP
		for j = 1:size(wrf(1).QCLOUD, 2)
			lwp(i, j) = trapz( squeeze(wrf(n).height( i, j, :, t )), squeeze(wrf(n).QCLOUD( i, j, :, t )) ) .* densityAir;
		end
	end
	h = pcolor( squeeze(wrf(1).XLONG(:,:,1)), squeeze(wrf(1).XLAT(:,:,1)), lwp); % Plot
	hold on
	plot(wrf(1).XLONG(25,25,1), wrf(1).XLAT(25,25,1), 'kx', 'markersize', 16); % Make X at center of domain
	% Border plot
	plot(wrf(1).XLONG(x,y(1),1), wrf(1).XLAT(x,y(1),1), 'k--'); plot(wrf(1).XLONG(x,y(end),1), wrf(1).XLAT(x,y(end),1), 'k--'); plot(wrf(1).XLONG(x(1),y,1), wrf(1).XLAT(x(1),y,1), 'k--'); plot(wrf(1).XLONG(x(end),y,1), wrf(1).XLAT(x(end),y,1), 'k--'); 
	% 				hold on
	% 				interval = 4:4:50;
	% 				quiver( squeeze(wrf(1).XLONG(interval, interval,1)), squeeze(wrf(1).XLAT(interval,interval,1)), squeeze(wrf(n).U(interval,interval, 10, t)), squeeze(wrf(n).V(interval,interval, 10, t)), 'r');
	% 		hold off
	set(h, 'edgecolor', 'none')
	% 		USBoundaries
	caxis([0 0.3])
	% 		xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
	colorbar
	title([name{1} ' LWP [kg/m^2]' 't = ' num2str( (t-1)./4 ) ' hr' ])
	
	% Right top
	axes('position', [0.78 0.655 0.2 sqPlotHeight])
	plot(timeVec, cv.avg.dthl)
	
	xlim([timeMarks(1), timeMarks(end)])
	set(gca,'xtick',timeMarks, 'xticklabel', [])
	grid on
	title('<\Delta\theta_l [K]>')
	% Right mid
	axes('position', [0.78 0.31 0.2 sqPlotHeight])
	plot(timeVec, cv.avg.dqt.*1000)
	
	xlim([timeMarks(1), timeMarks(end)])
	set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
	grid on
	title('<\Deltaq_t [g/kg]>')
	
	% Bottom 1
	axes('position', [0.04 0.04 longPlotWidth 0.2])
	plot(timeVec, cv.avg.shf, 'r--')
	hold on; plot(timeVec, cv.avg.lhf, 'b--')
	plot(les.timeTS./3600, les.sh, 'r')
	plot(les.timeTS./3600, les.lh, 'b')
	
	xlim([timeMarks(1), timeMarks(end)])
	set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
	grid on
	title('<Sfc fluxes [W/m^2]>')
	% Bottom 2
	axes('position', [0.04+longPlotWidth+offset 0.04 longPlotWidth 0.2])
	plot(timeVec, cv.avg.meanThetaL)
	hold on; plot(les.time./3600, les.meanTL, 'k')
	
	xlim([timeMarks(1), timeMarks(end)])
	set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
	grid on
	title('<\theta_{L,BL}>')
	% Bottom 3
	axes('position', [0.04+longPlotWidth*2+offset*2 0.04 longPlotWidth 0.2])
	plot(timeVec, cv.avg.meanQT.*1000)
	hold on; plot(les.time./3600, les.meanQT, 'k')
	
	xlim([timeMarks(1), timeMarks(end)])
	set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
	grid on
	title('<q_{t,BL}>')
	
	hgexport(gcf, [plotDir 'mapLWP_' sprintf('%.3d', t) name{1} '.png'], png)
	hgexport(gcf, [plotDir 'mapLWP_' sprintf('%.3d', t) name{1} '.eps'], eps)
	pause(0.1)
end

%%
figure
for t = 1:8:length(cv.time)
	h = pcolor( squeeze(wrf(1).XLONG(x,y,1)), squeeze(wrf(1).XLAT(x,y,1)), squeeze(cv.spatial.blRH( :, :, t ))); % Plot
	hold on
	plot(wrf(1).XLONG(25,25,1), wrf(1).XLAT(25,25,1), 'kx', 'markersize', 16); % Make X at center of domain
	% Border plot
	plot(wrf(1).XLONG(x,y(1),1), wrf(1).XLAT(x,y(1),1), 'k--'); plot(wrf(1).XLONG(x,y(end),1), wrf(1).XLAT(x,y(end),1), 'k--'); plot(wrf(1).XLONG(x(1),y,1), wrf(1).XLAT(x(1),y,1), 'k--'); plot(wrf(1).XLONG(x(end),y,1), wrf(1).XLAT(x(end),y,1), 'k--'); 
	% 				hold on
	% 				interval = 4:4:50;
	% 				quiver( squeeze(wrf(1).XLONG(interval, interval,1)), squeeze(wrf(1).XLAT(interval,interval,1)), squeeze(wrf(n).U(interval,interval, 10, t)), squeeze(wrf(n).V(interval,interval, 10, t)), 'r');
	% 		hold off
	set(h, 'edgecolor', 'none')
	% 		USBoundaries
	caxis([0.8 1])
			xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
	colorbar
	title([name{1} ' RH [-] ' 't = ' num2str( (t-1)./4 ) ' hr' ])
	hgexport(gcf, [plotDir 'mapRH_' sprintf('%.3d', t) name{1} '.png'], png)
	hgexport(gcf, [plotDir 'mapRH_' sprintf('%.3d', t) name{1} '.eps'], eps)
	pause(0.1)
end

return

%% CV Plot 2
i = 16; j = 16; % Temporary, but choose center of domain

sqPlotHeight = 0.3;
longPlotWidth = 0.28;
offset = 0.05;

timeMarks = ((1:16:length(cv.time)) - 1) ./ 4;
timeVec = (1:length(cv.time))./4;

figure('position', [120 20 1500 900])

% Left top
axes('position', [0.04 0.655 0.2 sqPlotHeight])
plot( timeVec, squeeze( cv.budgetHeat.w_e( i, j, : ) ), 'g' )
hold on, plot( timeVec, squeeze( cv.budgetMoisture.backedOutw_e( i, j, :) ), 'b' )
plot( timeVec, squeeze( cv.budgetHeat.backedOutw_e( i, j, :) ), 'r' )
plot( timeVec, squeeze( cv.budgetHeat.meanw_e(i, j, :) ), 'k')

xlim([timeMarks(1), timeMarks(end)])
ylim([-0.02, 0.02])
set(gca,'xtick',timeMarks, 'xticklabel', [])
grid on
title('w_e [m/s]')
legend('mass', 'MB', 'HB', 'mean', 'orientation', 'horizontal', 'location', 'best')

% Left mid
axes('position', [0.04 0.31 0.2 sqPlotHeight])
plot( timeVec, squeeze( cv.budgetHeat.wSub(i, j, :) ), 'm'); hold on
plot( timeVec, squeeze( cv.budgetHeat.wSubAbove(i, j, :) ), 'r--');


legend('w_s', 'w_{s+2pt}', 'location', 'best')

xlim([timeMarks(1), timeMarks(end)])
set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
grid on
title('w_{sub}')

% Center LWP
axes('position', [0.28 0.31 0.45 0.645])
clear qtMat
gridZ = 1:10:1200;
qlMat = zeros(length(gridZ), length(cv.time));

for t = 1:length(cv.time)
	idx = zeros(1, size(cv.heightAGL, 3));
	ql = zeros(1, size(cv.heightAGL, 3));
	for z = 1:size(cv.heightAGL, 3)
		idx(z) = findMatch(gridZ, cv.heightAGL( i, j, z, t ));
		% 		qlMat(idx(z),t) = cv.spatial.ql( i, j, z, t );
		ql(z) = cv.spatial.ql( i, j, z, t );
	end
	a = gridZ(idx); b = ql;
	a( ql == 0 ) = []; b( ql == 0 ) = [];
	if length(a) > 1
		qlInterp = interp1( a, b, gridZ );
		qlMat(:, t) = qlInterp;
	else
		qlMat(:, t) = 0;
	end
end

%%%%% NOTE TO SELF: COMBINE WITH ABOVE AND MAKE TWO MATS %%%%%

% % gridZ = linspace(1,1000);
% clear qtMat
% gridZ = 1:10:1000;
% qtMat = zeros(length(gridZ), length(cv.time));
%
% for t = 1:length(cv.time)
% 	idx = zeros(1, size(cv.heightAGL, 3));
% 	ql = zeros(1, size(cv.heightAGL, 3));
% 	for z = 1:size(cv.heightAGL, 3)
% 		idx(z) = findMatch(gridZ, cv.heightAGL( i, j, z, t ));
% 		% 		qlMat(idx(z),t) = cv.spatial.ql( i, j, z, t );
% 		ql(z) = cv.spatial.qt( i, j, z, t );
% 	end
% 	a = gridZ(idx); b = ql; c = diff(a) == 0;
% 	a( c ) = []; b( c ) = [];
% 	if length(a) > 1
% 		qlInterp = interp1( a, b, gridZ );
% 		qlMat(:, t) = qlInterp;
% 	else
% 		qlMat(:, t) = 0;
% 	end
% end
%%%%%

noCloud = squeeze(cv.spatial.zi(i, j, :)) == squeeze(cv.spatial.zb(i, j, :));
noCloudZi = squeeze(cv.spatial.zi(i, j, :)); noCloudZi( ~noCloud ) = NaN;
hold on
imagesc( timeVec, gridZ, qlMat )
plot( timeVec, squeeze(cv.spatial.zi(i, j, :)), 'w')
plot( timeVec, squeeze(cv.spatial.zb(i, j, :)), 'w--')
plot( timeVec, squeeze(cv.cloudBase.predictedzb(i, j, :)), 'g')
plot( timeVec, squeeze(cv.cloudBase.predictedzi(i, j, :)), 'g')
hold on; plot(les.timeTS./3600, les.zi1, 'k')
plot( les.timeTS./3600, les.zb, 'k--')
plot( timeVec, noCloudZi, 'r')

xlim([timeMarks(1), timeMarks(end)])
ylim([0 1200])
title([name{1} ' cloud boundaries'])
legend('WRF z_i', 'WRF z_b', 'MLM z_b', 'LES z_i', 'LES z_b', 'location', 'best')

% Right top
axes('position', [0.78 0.655 0.2 sqPlotHeight])
plot(timeVec, cv.avg.dthl)

xlim([timeMarks(1), timeMarks(end)])
set(gca,'xtick',timeMarks, 'xticklabel', [])
grid on
title('<\Delta\theta_l [K]>')
% Right mid
axes('position', [0.78 0.31 0.2 sqPlotHeight])
plot(timeVec, cv.avg.dqt.*1000)

xlim([timeMarks(1), timeMarks(end)])
set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
grid on
title('<\Deltaq_t [g/kg]>')

% Bottom 1
axes('position', [0.04 0.04 longPlotWidth 0.2])
hold on
plot(timeVec, squeeze(cv.cloudBase.dzb_dt_m1(i, j, :)), 'color', [0.7 0.7 0.7])
plot(timeVec, squeeze(cv.cloudBase.dzb_dt_m2(i, j, :)), 'g-*')
plot(timeVec, squeeze(cv.cloudBase.dzb_dt_m3(i, j, :)), 'y-o')
plot(timeVec, squeeze(cv.cloudBase.dzb_dt_m4(i, j, :)), 'b')
% plot(timeVec, squeeze(cv.cloudBase.dzb_dt_m5(i, j, :)), 'c-+')
axis tight
grid on
legend('advZi', 'advM', 'w_e', 'LHF', 'location', 'best', 'orientation', 'horizontal')

xlim([timeMarks(1), timeMarks(end)])
set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
grid on
title('Moisture budget')
ylabel('dz_b/dt [m/s]')
% Bottom 2
axes('position', [0.04+longPlotWidth+offset 0.04 longPlotWidth 0.2])
hold on
plot(timeVec, squeeze(cv.cloudBase.dzb_dt_h1(i, j, :)), 'color', [0.7 0.7 0.7])
plot(timeVec, squeeze(cv.cloudBase.dzb_dt_h2(i, j, :)), 'g-^')
plot(timeVec, squeeze(cv.cloudBase.dzb_dt_h3(i, j, :)), 'y-s')
plot(timeVec, squeeze(cv.cloudBase.dzb_dt_h4(i, j, :)), 'r')
plot(timeVec, squeeze(cv.cloudBase.dzb_dt_h5(i, j, :)), 'm--')
% plot(cv.time, squeeze(cv.cloudBase.dzb_dt_h6(i, j, :)), 'c-+')

axis tight
grid on
legend('advZi', 'advH', 'w_e', 'SHF', 'rad', 'location', 'best', 'orientation', 'horizontal')

xlim([timeMarks(1), timeMarks(end)])
set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
grid on
title('Heat budget')
% Bottom 3
axes('position', [0.04+longPlotWidth*2+offset*2 0.04 longPlotWidth 0.2])
plot( timeVec, squeeze( cv.cloudBase.dzb_dt(i, j, :) ), 'k' )
% hold on, plot( cv.time, squeeze( cv.cloudBase.dzb_dt_FD(i, j, :) ), 'g' )
hold on, plot( timeVec, squeeze( cv.cloudBase.dzb_dt_heat(i, j, :) ), '--r' )
plot( timeVec, squeeze( cv.cloudBase.dzb_dt_moisture(i, j, :) ), '--b' )
plot(timeVec, 0.2.*squeeze( cv.spatial.lwp(i, j, :) ), 'b^')

legend('ML theory', 'heat', 'moisture', 'location', 'best')
xlim([timeMarks(1), timeMarks(end)])
set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
grid on
title('dz_b/dt')
	
hgexport(gcf, [plotDir 'budget_' name{1} '.png'], png)
hgexport(gcf, [plotDir 'budget_' name{1} '.eps'], eps)

%% CV Plot 3 BL AVG
i = 16; j = 16; % Temporary, but choose center of domain

sqPlotHeight = 0.3;
longPlotWidth = 0.28;
offset = 0.05;

timeMarks = ((1:16:length(cv.time)) - 1) ./ 4;
timeVec = (1:length(cv.time))./4;

figure('position', [120 20 1500 900])

% Left top
axes('position', [0.04 0.655 0.2 sqPlotHeight])
plot( timeVec, squeeze( cv.avg.budgetHeat.w_e(: ) ), 'g' )
hold on, plot( timeVec, squeeze( cv.avg.budgetMoisture.backedOutw_e(:) ), 'b' )
plot( timeVec, squeeze( cv.avg.budgetHeat.backedOutw_e(:) ), 'r' )
plot( timeVec, squeeze( cv.avg.budgetHeat.meanw_e(:) ), 'k')

xlim([timeMarks(1), timeMarks(end)])
ylim([-0.02, 0.02])
set(gca,'xtick',timeMarks, 'xticklabel', [])
grid on
title('w_e [m/s]')
legend('mass', 'MB', 'HB', 'mean', 'orientation', 'horizontal', 'location', 'best')

% Left mid
axes('position', [0.04 0.31 0.2 sqPlotHeight])
plot( timeVec, squeeze( cv.avg.budgetHeat.wSub(:) ), 'm'); hold on
plot( timeVec, squeeze( cv.avg.budgetHeat.wSubAbove(:) ), 'r--');


legend('w_s', 'w_{s+2pt}', 'location', 'best')

xlim([timeMarks(1), timeMarks(end)])
set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
grid on
title('w_{sub}')

% Center LWP
axes('position', [0.28 0.31 0.45 0.645])
clear qtMat
gridZ = 1:10:1200;
qlMat = zeros(length(gridZ), length(cv.time));

for t = 1:length(cv.time)
	idx = zeros(1, size(cv.heightAGL, 3));
	ql = zeros(1, size(cv.heightAGL, 3));
	for z = 1:size(cv.heightAGL, 3)
		idx(z) = findMatch(gridZ, cv.heightAGL( i, j, z, t ));
		% 		qlMat(idx(z),t) = cv.spatial.ql( i, j, z, t );
		ql(z) = cv.avg.ql( z, t );
	end
	a = gridZ(idx); b = ql;
	a( ql == 0 ) = []; b( ql == 0 ) = [];
	if length(a) > 1
		qlInterp = interp1( a, b, gridZ );
		qlMat(:, t) = qlInterp;
	else
		qlMat(:, t) = 0;
	end
end

%%%%% NOTE TO SELF: COMBINE WITH ABOVE AND MAKE TWO MATS %%%%%

% % gridZ = linspace(1,1000);
% clear qtMat
% gridZ = 1:10:1000;
% qtMat = zeros(length(gridZ), length(cv.time));
%
% for t = 1:length(cv.time)
% 	idx = zeros(1, size(cv.heightAGL, 3));
% 	ql = zeros(1, size(cv.heightAGL, 3));
% 	for z = 1:size(cv.heightAGL, 3)
% 		idx(z) = findMatch(gridZ, cv.heightAGL( i, j, z, t ));
% 		% 		qlMat(idx(z),t) = cv.spatial.ql( i, j, z, t );
% 		ql(z) = cv.spatial.qt( i, j, z, t );
% 	end
% 	a = gridZ(idx); b = ql; c = diff(a) == 0;
% 	a( c ) = []; b( c ) = [];
% 	if length(a) > 1
% 		qlInterp = interp1( a, b, gridZ );
% 		qlMat(:, t) = qlInterp;
% 	else
% 		qlMat(:, t) = 0;
% 	end
% end
%%%%%

noCloud = squeeze(cv.avg.zi(:)) == squeeze(cv.avg.zb(:));
noCloudZi = squeeze(cv.avg.zi(:)); noCloudZi( ~noCloud ) = NaN;
hold on
imagesc( timeVec, gridZ, qlMat )
plot( timeVec, squeeze(cv.avg.zi(:)), 'w')
plot( timeVec, squeeze(cv.avg.zb(:)), 'w--')
plot( timeVec, squeeze(cv.avg.cloudBase.predictedzb(:)), 'g')
plot( timeVec, squeeze(cv.avg.cloudBase.predictedzi(:)), 'g')
hold on; plot(les.timeTS./3600, les.zi1, 'k')
plot( les.timeTS./3600, les.zb, 'k--')
plot( timeVec, noCloudZi, 'r')

xlim([timeMarks(1), timeMarks(end)])
ylim([0 1200])
title([name{1} ' cloud boundaries'])
legend('WRF z_i', 'WRF z_b', 'MLM z_b', 'LES z_i', 'LES z_b', 'location', 'best')

% Right top
axes('position', [0.78 0.655 0.2 sqPlotHeight])
plot(timeVec, cv.avg.dthl)

xlim([timeMarks(1), timeMarks(end)])
set(gca,'xtick',timeMarks, 'xticklabel', [])
grid on
title('<\Delta\theta_l [K]>')
% Right mid
axes('position', [0.78 0.31 0.2 sqPlotHeight])
plot(timeVec, cv.avg.dqt.*1000)

xlim([timeMarks(1), timeMarks(end)])
set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
grid on
title('<\Deltaq_t [g/kg]>')

% Bottom 1
axes('position', [0.04 0.04 longPlotWidth 0.2])
hold on
plot(timeVec, squeeze(cv.avg.cloudBase.dzb_dt_m1(:)), 'color', [0.7 0.7 0.7])
plot(timeVec, squeeze(cv.avg.cloudBase.dzb_dt_m2(:)), 'g-*')
plot(timeVec, squeeze(cv.avg.cloudBase.dzb_dt_m3(:)), 'y-o')
plot(timeVec, squeeze(cv.avg.cloudBase.dzb_dt_m4(:)), 'b')
% plot(timeVec, squeeze(cv.cloudBase.dzb_dt_m5(i, j, :)), 'c-+')
axis tight
grid on
legend('advZi', 'advM', 'w_e', 'LHF', 'location', 'best', 'orientation', 'horizontal')

xlim([timeMarks(1), timeMarks(end)])
set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
grid on
title('Moisture budget')
ylabel('dz_b/dt [m/s]')
% Bottom 2
axes('position', [0.04+longPlotWidth+offset 0.04 longPlotWidth 0.2])
hold on
plot(timeVec, squeeze(cv.avg.cloudBase.dzb_dt_h1(:)), 'color', [0.7 0.7 0.7])
plot(timeVec, squeeze(cv.avg.cloudBase.dzb_dt_h2(:)), 'g-^')
plot(timeVec, squeeze(cv.avg.cloudBase.dzb_dt_h3(:)), 'y-s')
plot(timeVec, squeeze(cv.avg.cloudBase.dzb_dt_h4(:)), 'r')
plot(timeVec, squeeze(cv.avg.cloudBase.dzb_dt_h5(:)), 'm--')
% plot(cv.time, squeeze(cv.cloudBase.dzb_dt_h6(i, j, :)), 'c-+')

axis tight
grid on
legend('advZi', 'advH', 'w_e', 'SHF', 'rad', 'location', 'best', 'orientation', 'horizontal')

xlim([timeMarks(1), timeMarks(end)])
set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
grid on
title('Heat budget')
% Bottom 3
axes('position', [0.04+longPlotWidth*2+offset*2 0.04 longPlotWidth 0.2])
plot( timeVec, squeeze( cv.avg.cloudBase.dzb_dt(:) ), 'k' )
% hold on, plot( cv.time, squeeze( cv.avg.cloudBase.dzb_dt_FD(:) ), 'g' )
hold on, plot( timeVec, squeeze( cv.avg.cloudBase.dzb_dt_heat(:) ), '--r' )
plot( timeVec, squeeze( cv.avg.cloudBase.dzb_dt_moisture(:) ), '--b' )
plot(timeVec, 0.2.*squeeze( cv.avg.lwp(:) ), 'b^')

legend('ML theory', 'heat', 'moisture', 'location', 'best')
xlim([timeMarks(1), timeMarks(end)])
set(gca,'xtick',timeMarks, 'xticklabel', timeMarks)
grid on
title('dz_b/dt')
	
hgexport(gcf, [plotDir 'budget_avg_' name{1} '.png'], png)
hgexport(gcf, [plotDir 'budget_avg_' name{1} '.eps'], eps)

%% TO DO LIST
% Make profile evolution plot
% Compute wsub using model
% - Compute dzi/dt using wsub
% - Get MLM zi
% Compute advection and budget terms using CV
% Convert this script into a function and run spread!
	
% %% BREAK
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	
% %% BREAK
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Load LWP
% % Domain averaged? Take single point for now
% i = 25; j = 25;
% 
% %%
% n1 = 1; n2 = 2; % Compare between configurations
% i = 25; j = 25;
% 
% for n = 1:length(wrf)
% 	wrf(n).QT = wrf(n).QCLOUD + wrf(n).QVAPOR;% + wrf(n).QRAIN;
% end
% 
% % plotVars = {'TL', 'QCLOUD', 'SH3D', 'EL_PBL'}; heightVars = {'height', 'height', 'height', 'heightstag'}; x = [285, 310; 0, 0.001; 0, 1.2; 0 500];
% % plotVars = {'U', 'V', 'EXCH_H', 'QCLOUD'};
% % plotVars = {'T', 'QT', 'QCLOUD', 'U'}; heightVars = {'height', 'height', 'height', 'height'};
% % % % % % % % % % % % % % % % plotVars = {'EL_PBL', 'QKE', 'SH3D', 'EXCH_H'};
% % % % % % % % % % % % % % % % heightVars = {'heightstag', 'height', 'height', 'heightstag'};
% 
% % plotVars = {'T', 'QCLOUD', 'QVAPOR', 'EL_PBL'}; heightVars = {'height', 'height', 'height', 'heightstag'}; x = [285, 310; 0, 0.001; 0, 0.1; 0 500];
% % plotVars = {'TL', 'QCLOUD', 'RTHRATEN', 'QKE'}; heightVars = {'height', 'height', 'height', 'height'}; x = [285, 310; 0, 0.001; -0.0005, 0.0005; 0 1.5];
% plotVars = {'TL', 'QCLOUD', 'EXCH_H', 'QKE'}; heightVars = {'height', 'height', 'heightstag', 'height'}; x = [285, 310; 0, 0.001; 0, 75; 0 1];
% % plotVars = {'TL', 'QCLOUD', 'EXCH_H', 'EL_PBL'}; heightVars = {'height', 'height', 'heightstag', 'heightstag'}; x = [285, 310; 0, 0.001; 0, 300; 0 500];
% 
% plotVars = {'TL', 'QCLOUD', 'EXCH_H', 'W'}; heightVars = {'height', 'height', 'heightstag', 'heightstag'}; x = [285, 310; 0, 0.001; 0, 75; -0.005 0.005];
% 
% yTop = 3000;
% % x = [280, 305; 0, 3000; 0, 2; 0 1500];
% % x = [-15, 15; -15, 15; 0, 300; 0 .005];
% % x = [285, 310; 0, 0.015; 0, .005; -5 8];
% % % % % % % % % % % % % % % x = [0, 200; 0, 1; 0, 1.5; 0 100];
% 
% figure('position', [50 50 1400 900])
% % for t = 1:480:length(wrf(1).time)
% for t = 1:4:length(wrf(1).time)
% 	for v = 1:length(plotVars)
% 		
% 		% Top
% 		subplot(2,length(plotVars),v)
% 		plot( squeeze(wrf(n1).(plotVars{v})(i, j, :, t)), squeeze(wrf(n1).(heightVars{v})(i, j, :, t)) );
% 		hold on
% 		plot( [-1000 1000], [wrf(n1).PBLH(i,j,t), wrf(n1).PBLH(i,j,t)], '--k')
% 		hold off
% 		title(plotVars{v})
% 		ylim([0 yTop])
% 		xlim(x(v,:))
% 		grid on
% 		
% 		% Bottom
% 		subplot(2,length(plotVars),length(plotVars) + v)
% 		plot( squeeze(wrf(n2).(plotVars{v})(i, j, :, t)), squeeze(wrf(n2).(heightVars{v})(i, j, :, t)) );
% 		hold on
% 		plot( [-1000 1000], [wrf(n2).PBLH(i,j,t), wrf(n2).PBLH(i,j,t)], '--k')
% 		hold off
% 		title(plotVars{v})
% 		ylim([0 yTop])
% 		xlim(x(v,:))
% 		grid on
% 		
% 	end
% 	text(0.4 .* x(4,2), 0.95.*yTop, ['t = ' num2str(t)])
% 	text(0.2 .* x(4,2), 0.8.*yTop, name{n1})
% 	text(0.2 .* x(4,2), 0.6.*yTop, name{n2})
% 	text(0.2 .* x(4,2), 0.4.*yTop, ['lwp1=' num2str(wrf(n1).lwp(i,j,t))])
% 	text(0.2 .* x(4,2), 0.3.*yTop, ['lwp2=' num2str(wrf(n2).lwp(i,j,t))])
% 	hgexport(gcf, ['~/Dropbox/tempPlots/' 'THERMOprofiles_' sprintf('%.2d', t) name{1} '_' name{2} '.png'], png)
% 	pause()
% end
% 
% %% LWP Map
% 
% wrf(1).XLONG = wrfLoad(path{1}, {'XLONG'});
% wrf(1).XLAT = wrfLoad(path{1}, {'XLAT'});
% 
% figure('position', [50 50 1400 900])
% for t = 1:4:size(wrf(n).T, 4)
% 	for n = 1:length(wrf)
% 		for i = 1:size(wrf(1).QCLOUD, 1)
% 			for j = 1:size(wrf(1).QCLOUD, 2)
% 				lwp(i, j) = trapz( squeeze(wrf(n).height( i, j, :, t )), squeeze(wrf(n).QCLOUD( i, j, :, t )) ) .* densityAir;
% 			end
% 		end
% 		subplot(1,2,n)
% 		
% 		h = pcolor( squeeze(wrf(1).XLONG(:,:,1)), squeeze(wrf(1).XLAT(:,:,1)), lwp);
% 		hold on
% 		interval = 4:4:50;
% 		quiver( squeeze(wrf(1).XLONG(interval, interval,1)), squeeze(wrf(1).XLAT(interval,interval,1)), squeeze(wrf(n).U(interval,interval, 10, t)), squeeze(wrf(n).V(interval,interval, 10, t)), 'r');
% 		hold off
% 		set(h, 'edgecolor', 'none')
% % 		USBoundaries
% 		caxis([0 0.3])
% 		xlabel('Longitude [deg]'), ylabel('Latitude [deg]'), title('LWP')
% 	end
% 	title(datestr(wrf(1).time(t)))
% 	colorbar
% 	hgexport(gcf, ['~/Dropbox/tempPlots/' 'mapLWP_' sprintf('%.3d', t) name{1} '_' name{2} '.png'], png)
% 	pause(0.1)
% 	clf
% end
% 
% return
% 
% % hgexport(gcf, ['tempPlots/' 'LWP_' name{1} '_' name{2} '.png'], png)
% 
% %% TKE BUDGET
% n1 = 1; n2 = 2; % Compare between configurations
% i = 25; j = 25;
% 
% yTop = 2000;
% 
% figure('position', [50 50 1400 900])
% % for t = 1:480:length(wrf(1).time)
% for t = 1:1:length(wrf(1).time)
% 	for v = 1:length(plotVars)
% 		% Top
% 		clf
% 		hold on
% 		plot(squeeze(wrf(n1).QSHEAR(i, j, :, t)), squeeze(wrf(n1).heightstag(i, j, :, t)), '--g')
% 		plot(squeeze(wrf(n1).QBUOY(i, j, :, t)), squeeze(wrf(n1).heightstag(i, j, :, t)), 'y--')
% 		plot(squeeze(wrf(n1).QDISS(i, j, :, t)), squeeze(wrf(n1).heightstag(i, j, :, t)), 'r--')
% 		plot(squeeze(wrf(n1).QWT(i, j, :, t)), squeeze(wrf(n1).heightstag(i, j, :, t)), 'b--')
% 		plot(squeeze(wrf(n1).DTKE(i, j, :, t)), squeeze(wrf(n1).height(i, j, :, t)), 'k--')
% 		hold off
% 		
% 		ylim([0 yTop])
% 		title(num2str(t))
% 		grid on
% 		
% 	end
% 	% 	hgexport(gcf, ['~/Dropbox/Presentations/SCM-VertAdv/' '_profile_' sprintf('%.2d', t) name1 '_' name2 '.png'], png)
% 	pause()
% end
% 
% %% Four panel TKE budget
% n1 = 1; n2 = 2;
% hoursEnd = 1:4;
% yTop = 1200;
% x1 = -0.03; x2 = 0.03;
% 
% figure('position', [50 50 1400 900])
% for h = 1:length(hoursEnd)
% 	
% 	resolution = 15; % Minutes
% 	endIdx = 1 + 60 ./ resolution .* hoursEnd(h); startIdx = endIdx - 60 ./ resolution; 
% 	
% 	subplot(2,length(hoursEnd),h)
% 	plot( mean( squeeze( wrf(n1).QBUOY(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n1).heightstag(i, j, :, startIdx:endIdx) ), 2), 'color', [1 0.6 0])
% 	hold on
% 	plot( mean( squeeze( wrf(n1).QSHEAR(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n1).heightstag(i, j, :, startIdx:endIdx) ), 2), 'g--')
% 	plot( mean( squeeze( wrf(n1).QDISS(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n1).heightstag(i, j, :, startIdx:endIdx) ), 2), 'r--')
% 	plot( mean( squeeze( wrf(n1).QWT(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n1).heightstag(i, j, :, startIdx:endIdx) ), 2), 'b--')
% 	
% 	ylim([0 yTop])
% 	xlim([x1 x2])
% 	xlabel('m2 s-2')
% 	
% 	subplot(2,length(hoursEnd),length(hoursEnd) + h)
% 	plot( mean( squeeze( wrf(n2).QBUOY(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n2).heightstag(i, j, :, startIdx:endIdx) ), 2), 'color', [1 0.6 0])
% 	hold on
% 	plot( mean( squeeze( wrf(n2).QSHEAR(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n2).heightstag(i, j, :, startIdx:endIdx) ), 2), 'g--')
% 	plot( mean( squeeze( wrf(n2).QDISS(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n2).heightstag(i, j, :, startIdx:endIdx) ), 2), 'r--')
% 	plot( mean( squeeze( wrf(n2).QWT(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n2).heightstag(i, j, :, startIdx:endIdx) ), 2), 'b--')
% 	title(['t = ' num2str(h) ' hr'])
% 	ylim([0 yTop])
% 	xlim([x1 x2])
% 	
% 	hgexport(gcf, ['~/Dropbox/tempPlots/' 'TKE_' sprintf('%.2d', t) name{1} '_' name{2} '.png'], png)
% 	
% end
% 
% 
% %% Normalize and time average
% 
% zRange = 0:0.05:2;
% 
% n = 1; % n of MYNN case
% i = 25; j = 25;
% 
% zinv = 800;
% 
% for t = 2:size(wrf(1).T, 4)
% 	haveDTKE = squeeze(wrf(n).DTKE(i, j, :, t));
% 	haveQKE = squeeze(wrf(n).QKE(i, j, :, t));
% 	haveEL = squeeze(wrf(n).EL_PBL(i, j, :, t));
% 	haveS = squeeze(wrf(n).QSHEAR(i, j, :, t)) ./ haveEL;
% 	haveB = squeeze(wrf(n).QBUOY(i, j, :, t)) ./ haveEL;
% 	haveD = squeeze(wrf(n).QDISS(i, j, :, t)) ./ haveEL;
% 	haveT = squeeze(wrf(n).QWT(i, j, :, t)) ./ haveEL;
% 	
% % 	haveTSQ = squeeze(wrf(n).TSQ(i, j, :, t));
% % 	haveQSQ = squeeze(wrf(n).QSQ(i, j, :, t));
% % 	haveCOV = squeeze(wrf(n).COV(i, j, :, t));
% % 	haveWSQ = squeeze(wrf(n).WSQ(i, j, :, t));
% 	
% 	haveSH3D = squeeze(wrf(n).SH3D(i, j, :, t));
% 	haveEXCH_H = squeeze(wrf(n).EXCH_H(i, j, :, t));
% 	
% % 	haveZ = squeeze(wrf(n).height(i,j,:,t) ./ wrf(n).PBLH(i, j, t));
% % 	haveZstag = squeeze(wrf(n).heightstag(i,j,:,t) ./ wrf(n).PBLH(i, j, t));
% 
% 	haveZ = squeeze(wrf(n).height(i,j,:,t) ./ zinv);
% 	haveZstag = squeeze(wrf(n).heightstag(i,j,:,t) ./ zinv);
% 	
% 	wrf(n).iDTKE(i, j, :, t) = interp1(haveZ, haveDTKE, zRange);
% 	wrf(n).iQKE(i, j, :, t) = interp1(haveZ, haveQKE, zRange);
% 	wrf(n).iS(i, j, :, t) = interp1(haveZstag, haveS, zRange);
% 	wrf(n).iB(i, j, :, t) = interp1(haveZstag, haveB, zRange);
% 	wrf(n).iD(i, j, :, t) = interp1(haveZstag, haveD, zRange);
% 	wrf(n).iT(i, j, :, t) = interp1(haveZstag, haveT, zRange);
% 	wrf(n).iEL(i, j, :, t) = interp1(haveZstag, haveEL, zRange);
% 	
% % 	wrf(n).iTSQ(i, j, :, t) = interp1(haveZ, haveTSQ, zRange);
% % 	wrf(n).iQSQ(i, j, :, t) = interp1(haveZ, haveQSQ, zRange);
% % 	wrf(n).iCOV(i, j, :, t) = interp1(haveZ, haveCOV, zRange);
% % 	wrf(n).iWSQ(i, j, :, t) = interp1(haveZ, haveWSQ, zRange);
% 	
% 	wrf(n).iSH3D(i, j, :, t) = interp1(haveZ, haveSH3D, zRange);
% 	wrf(n).iEXCH_H(i, j, :, t) = interp1(haveZstag, haveEXCH_H, zRange);
% 	
% end
% 
% 
% % 2 - 4
% tRange =  (1+(2*4)):(4*4 + 1); % 4-6
% lesRange = 2:4;
% % 
% % 11 - 13
% % tRange =  (1+(11*60)):(13*60 + 1); % 11-13
% % lesRange = 11:13;
% 
% 
% for z = 1:length(zRange)
% 	avgDTKE(n,z) = mean(wrf(n).iDTKE(i,j,z,tRange), 4);
% 	avgQKE(n,z) = mean(wrf(n).iQKE(i,j,z,tRange), 4);
% 	avgS(n,z) = mean(wrf(n).iS(i,j,z,tRange), 4);
% 	avgB(n,z) = mean(wrf(n).iB(i,j,z,tRange), 4);
% 	avgD(n,z) = mean(wrf(n).iD(i,j,z,tRange), 4);
% 	avgT(n,z) = mean(wrf(n).iT(i,j,z,tRange), 4);
% 	avgEL(n,z) = mean(wrf(n).iEL(i,j,z,tRange), 4);
% 	
% % 	avgTSQ(n,z) = mean(wrf(n).iTSQ(i,j,z,tRange), 4);
% % 	avgQSQ(n,z) = mean(wrf(n).iQSQ(i,j,z,tRange), 4);
% % 	avgCOV(n,z) = mean(wrf(n).iCOV(i,j,z,tRange), 4);
% % 	avgWSQ(n,z) = mean(wrf(n).iWSQ(i,j,z,tRange), 4);
% % 	avgWSQms(n,z) = mean(wrf(n).iWSQ(i,j,z,tRange) .* wrf(n).iQKE(i,j,z,tRange), 4);
% 	
% 	avgSH3D(n,z) = mean(wrf(n).iSH3D(i,j,z,tRange), 4);
% 	avgEXCH_H(n,z) = mean(wrf(n).iEXCH_H(i,j,z,tRange), 4);
% 	
% end
% 
% avgZB(n) = mean(wrf(n).zb(i,j,tRange));
% 
% % Average LES
% startIdx = 1 + lesRange(1)*3600/20/60; endIdx = 1 + lesRange(end)*3600/20/60;
% startIdx2 = 1 + lesRange(1)*3600/20; endIdx2 = 1 + lesRange(end)*3600/20;
% midIdx = floor( (startIdx2 + endIdx2) / 2);
% for z = 1:size(les.wthl,1)
% 	les.avgTKE(z) = mean(les.tke(z,startIdx:endIdx));
% 	les.avgWTHL(z) = mean(les.wthl(z,startIdx:endIdx));
% 	les.avgWQT(z) = mean(les.wqt(z,startIdx:endIdx));
% 	les.avgWW(z) = mean(les.ww(z,startIdx:endIdx));
% 	les.avgBUOY(z) = mean(les.buoy(z,startIdx:endIdx));
% 	les.avgSHEAR(z) = mean(les.shr(z,startIdx:endIdx));
% 	les.avgTRANS(z) = mean(les.trans(z,startIdx:endIdx));
% 	les.avgDISS(z) = mean(les.diss(z,startIdx:endIdx));
% end
% les.avgZC = mean(les.zc(startIdx2:endIdx2));
% les.avgZB = mean(les.zb(startIdx2:endIdx2));
% norm = 1 ./ les.zi1(floor( (startIdx2 + endIdx2) / 2 ));
% 
% 
% %% Plot
% %%%%%%%%%%% 2-hr avg %%%%%%%%%%%
% figure('position', [50 50 1200 600])
% 
% subplot(1,3,1)
% plot( .5 .* avgQKE(n,:), zRange, 'k' ) % 2*TKE
% hold on
% plot(les.avgTKE(2:end), les.zt(2:end) .* norm, 'k--')
% 
% plot([-1000 1000], [1, 1], 'k:', 'linewidth', 1) % PBLH
% plot([-1000 1000], [les.avgZB.*norm, les.avgZB.*norm], 'k:', 'linewidth', 1) % PBLH
% plot([-1000 1000], [avgZB(1)./zinv, avgZB(1)./zinv], 'r-.', 'linewidth', 1) % zb MYNN
% plot([0 0], [0 10000], 'k:') % 0 origin line
% legend('TKE', 'TKE_{les}')
% 
% xlim([0 1])
% ylim([0 1.5])
% xlabel('TKE [m2 s-2]')
% ylabel('z / z_i [-]')
% title TKE
% 
% subplot(1,3,2)
% % plot( 100 .* avgDTKE(n,:), zRange, 'k' ) % DTKE
% hold on
% % plot( 100.* squeeze(wrf(n).TKE_PBL(i, j, :, t)), squeeze(wrf(n).heightstag(i, j, :, t)), 'k^--' ) % TKE FROM PBL -- VERIFIED MATCH
% plot( 100.* avgS(n,:), zRange, 'g' ) % Shear production
% plot( 100.* avgB(n,:), zRange, 'm' ) % Buoyant production
% plot(100.*les.avgSHEAR(1:2:end), les.zt(1:2:end) .* norm, 'go', 'markersize', 2)
% plot(100.*les.avgBUOY(2:2:end), les.zm(2:2:end) .* norm, 'mo', 'markersize', 2)
% 
% % divide by L
% 
% plot([-1000 1000], [1, 1], 'k:', 'linewidth', 1) % PBLH
% % plot([-1000 1000], [les.avgZC.*norm, les.avgZC.*norm], 'k', 'linewidth', 1) % PBLH
% plot([-1000 1000], [les.avgZB.*norm, les.avgZB.*norm], 'k:', 'linewidth', 1) % zb LES
% plot([-1000 1000], [avgZB(1)./zinv, avgZB(1)./zinv], 'r-.', 'linewidth', 1) % zb MYNN
% plot([0 0], [0 10000], 'k:') % 0 origin line
% legend('S', 'B', 'S_{les}', 'B_{les}')
% 
% ylim([0 1.5])
% xlim([-0.1 0.45])
% xlabel('[cm2 s-3]')
% title Production
% 
% subplot(1,3,3)
% 
% plot(-100.* avgD(n,:), zRange, 'r' ) % Dissipation
% hold on
% plot(100.* avgT(n,:), zRange, 'b' ) % Transport
% plot(-100.*les.avgDISS(1:2:end), les.zm(1:2:end) .* norm, 'ro', 'markersize', 2)
% plot(100.*les.avgTRANS(1:2:end), les.zm(1:2:end) .* norm, 'bo', 'markersize', 2)
% 
% plot([-1000 1000], [1, 1], 'k:', 'linewidth', 1) % PBLH
% plot([-1000 1000], [les.avgZB.*norm, les.avgZB.*norm], 'k:', 'linewidth', 1) % PBLH
% plot([-1000 1000], [avgZB(1)./zinv, avgZB(1)./zinv], 'r-.', 'linewidth', 1) % zb MYNN
% plot([0 0], [0 10000], 'k:') % 0 origin line
% legend('D', 'T', 'D_{les}', 'T_{les}', 'location', 'southwest')
% 
% ylim([0 1.5])
% xlim([-0.7 0.4])
% xlabel('[cm2 s-3]')
% title ('Dissip. & Transport')
% 
% %% INT TKE
% 
% n = 1;
% i = 25; j = 25;
% for n = 1:length(wrf)
% 	for t = 1:length(wrf(n).time)
% 		% 	wrf(n).intTKE(1, 1, t) = sum( squeeze(wrf(n).TKE_PBL( 1, 1, 1:70, t)) .* squeeze( wrf(n).DNW( 1:70, t ) ) ) ./ sum( wrf(n).DNW( 1:70, t ) );
% 		wrf(n).intTKE(i, j, t) = trapz( squeeze(wrf(n).height( i, j, :, t )), 0.5 .* squeeze(wrf(n).QKE( i, j, :, t )) );
% 	end
% end
% % for t = 1:length(les.time)
% % 	les.TKEintManual(t) = trapz(les.zt,les.tke(:,t));
% % end
% 
% figure; 
% plot( (1:length(wrf(1).time)) ./ 4, squeeze(wrf(1).intTKE(i,j,:)), 'r'); hold on
% plot( (1:length(wrf(2).time)) ./ 4, squeeze(wrf(2).intTKE(i,j,:)), 'b');
% % plot( (1:length(wrf(3).time)) ./ 60, squeeze(wrf(3).intTKE(i,j,:)), 'g');
% hold on
% plot(les.timeTS./3600, les.tkeint, 'k')
% % plot(les.time./3600, les.TKEintManual, 'r--')
% % legend('MYNN', 'LES', 'location', 'best')
% h = legend(name{1}, name{2}, 'LES', 'location', 'best');
% legend('boxoff');
% set(gca, 'xtick', 0:2:24)
% xlim([0 24])
% grid on
% ylabel('Vert. int. TKE [m3 s-2]')
% xlabel('Time [hr]')
% 
% hgexport(gcf, ['~/Dropbox/tempPlots/' 'INTTKE_' name{1} '_' name{2} '.png'], png)
end