%% Compare WRF ideal simulations

clear all
clc
close all

%% Define physics constants
defineConstants;
%%

name{1} = 'MYNN-3.6';
path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.6/wrfout_d01_20150710_0800.nc';
% % name{2} = 'MYNN-3.6-NBL';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.6noBouLac/wrfout_d01_20150710_0800.nc';
% % name{2} = 'MYNN-3.6-MYNN3';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.6MYNN3/wrfout_d01_20150710_0800.nc';
% % % % % name{1} = 'MYNN-3.7.1';
% % % % % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1mynn/wrfout_d01_20150710_0800.nc';
% % % % % name{2} = 'MYNN-3.7.1-0F';
% % % % % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1_equator/wrfout_d01_20150710_0800.nc';

% name{1} = 'MYNN-3.7.1';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1mynn_0f/wrfout_d01_20150710_0800.nc';
% name{1} = 'MYNN-3.8.1-ML';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1mynnML2_0f/wrfout_d01_20150710_0800.nc';
% name{2} = 'MYNN-3.8.1-0f';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1mynn_0f/wrfout_d01_20150710_0800.nc';
% 
% 
% name{1} = 'YSU-TD-3.7.1';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1ysutd_0f/wrfout_d01_20150710_0800.nc';
% name{2} = 'YSU-BUOY-3.7.1';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1ysubuoy_0f/wrfout_d01_20150710_0800.nc';
% name{3} = 'YSU-TD-3.8.1';
% path{3} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1ysutd_0f/wrfout_d01_20150710_0800.nc';


% % % name{2} = 'MYNN-3.8.1';
% % % % % % % % % % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1NEWWIND/wrfout_d01_20150710_0800.nc';
% % % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1/wrfout_d01_20150710_0800.nc';
% name{1} = 'MYNN-3.8.1';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1_24hr/wrfout_d01_20150710_0800.nc';
% name{1} = 'MYNN-3.8.1-TEMF';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1temf/wrfout_d01_20150710_0800.nc';

% name{2} = 'YSU-3.6';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.6ysu/wrfout_d01_20150710_0800.nc';
% name{1} = 'YSU-TD-3.8.1';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1ysuTD/wrfout_d01_20150710_0800.nc';
% name{1} = 'YSU-TD-3.7.1';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1ysuTD/wrfout_d01_20150710_0800.nc';
% % name{2} = 'YSU-BY-3.7.1';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1ysubuoy/wrfout_d01_20150710_0800.nc';
% name{2} = 'YSU-3.8.1';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1ysu/wrfout_d01_20150710_0800.nc';

%% HY test
% name{1} = 'MYNN-3.9-GB';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01GB/wrfout_d01_20150710_0800.nc';
% 
name{1} = 'MYNN-3.9-GB-51';
path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfSCM/testGB_39vanillathl_51/wrfout_d01_2015-07-10_08:00:00';

% 
% % name{2} = 'MYNN-3.8.1';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1mynn_0f/wrfout_d01_20150710_0800.nc';
% name{2} = 'MYNN-3.8.1-51';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1mynn_51/wrfout_d01_20150710_0800.nc';
% 
% name{2} = 'MYNN-3.9';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn/wrfout_d01_20150710_0800.nc';


%%%%%%%%%%%%%%%

% % Rad 0
% name{1} = 'MYNN-3.9-GB-0r';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_0rad/wrfout_d01_2015-07-10_08_00_00';
% name{2} = 'MYNN-3.9-0r';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51_norad/wrfout_d01_20150710_0800.nc';
% 
% % Rad 5
name{1} = 'MYNN-3.9-GB-5r';
path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_5rad/wrfout_d01_2015-07-10_08_00_00';
name{2} = 'MYNN-3.9-GB-4r';
path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_4rad/wrfout_d01_2015-07-10_08_00_00';

%%
% % Rad 5
% name{1} = 'MYNN-3.9-GB-4r-1icloud';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_4rad/wrfout_d01_2015-07-10_08_00_00';
name{1} = 'MYNN-3.9-GB-4r-0icloud';
path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_4rad_0icloudbl/wrfout_d01_2015-07-10_08_00_00';

% ELB
% name{1} = 'MYNN-3.9-GB-200ELP';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/200ELPROFILE_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% name{1} = 'MYNN-3.9-GB-2ELB';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/2ELB_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% 
% name{1} = 'MYNN-3.9-GB-EDMF1';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
name{2} = 'MYNN-3.9-GB-10ELPROFILE';
path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/100ELPROFILE_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00'; % actually 50 EL profile

% 2x mixing (EL and K)
% name{1} = 'MYNN-3.9-GB-2EL';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/2EL_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% name{2} = 'MYNN-3.9-GB-0.5EL';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/0.5EL_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% name{1} = 'MYNN-3.9-GB-4r-0icloud';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_4rad_0icloudbl/wrfout_d01_2015-07-10_08_00_00';
% name{2} = 'MYNN-3.9-GB-2K';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/2K_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';


% name{2} = 'MYNN-3.9-GB-4r-1icloud-1ml';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_4rad_1icloudbl_1ml/wrfout_d01_2015-07-10_08_00_00';
% name{2} = 'MYNN-3.9-GB-5r-1icloud-1ml';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_5rad_1icloudbl_1ml/wrfout_d01_2015-07-10_08_00_00';

% name{2} = 'MYNN-3.9-5r';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51_5rad/wrfout_d01_2015-07-10_08_00_00';
% 
% % 3.9 no GB
% name{1} = 'MYNN-3.9-0r';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51_norad/wrfout_d01_20150710_0800.nc';
% name{2} = 'MYNN-3.9-5r';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51_5rad/wrfout_d01_2015-07-10_08_00_00';
% 
% name{1} = 'MYNN-3.9-GB-0r';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_0rad/wrfout_d01_2015-07-10_08_00_00';
% name{2} = 'MYNN-3.9-GB-5r';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_5rad/wrfout_d01_2015-07-10_08_00_00';

%% BLOCK DATA SET DYCOMS RF01 EXPERIMENT

dycomsrf01 = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/';

% % Control 3.6
% name{1} = 'MYNN-3.6';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.6/wrfout_d01_20150710_0800.nc';
% 
% % Control 3.7.1
% name{1} = 'MYNN-3.7.1';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1mynn_0f/wrfout_d01_20150710_0800.nc';

% CONTROL CASE NO EDMF
name{1} = 'MYNN-3.9';
path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_4rad_0icloudbl/wrfout_d01_2015-07-10_08_00_00'; % CONTROL
% 
% % CONTROL CASE WITH EDMF = 1
name{2} = 'MYNN-3.9-EDMF';
path{2} = [dycomsrf01 'edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% 
% % 0.5 EL
% name{1} = 'MYNN-3.9-0.5EL';
% path{1} = [dycomsrf01 '0.5EL_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % 0.5 EL PROFILE
% name{2} = 'MYNN-3.9-0.5ELPROF';
% path{2} = [dycomsrf01 '0.5ELPROFILE_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % 
% % 
% % 2.0 EL
% name{1} = 'MYNN-3.9-1.5EL';
% path{1} = [dycomsrf01 '1.5_EL_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % 2.0 EL PROFILE
% name{2} = 'MYNN-3.9-1.5ELPROF';
% path{2} = [dycomsrf01 '1.5_ELPROFILE_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % 
% % % 0.5 EL
% % 
% name{1} = 'MYNN-3.9-EDMF-0.5EL';
% path{1} = [dycomsrf01 '0.5EL_edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % 0.5 EL PROFILE
% name{2} = 'MYNN-3.9-EDMF-0.5ELPROF';
% path{2} = [dycomsrf01 '0.5ELPROFILE_edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % 
% % 
% % 2.0 EL
% name{1} = 'MYNN-3.9-EDMF-1.5EL';
% path{1} = [dycomsrf01 '1.5_EL_edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % 2.0 EL PROFILE
% name{2} = 'MYNN-3.9-EDMF-1.5ELPROF';
% path{2} = [dycomsrf01 '1.5_ELPROFILE_edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];

% NOBDY
% name{1} = 'YSU-3.9-TD-0BDY';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/nobdy_edmf0_pbl1_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';

% name{1} = 'MYNN-3.9-GB-0BDY';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/nobdy_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';

% WITH SUBS

% name{1} = 'YSU-3.9-TD-0BDY';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/375Dsubs_edmf0_pbl1_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';

% name{2} = 'MYNN-3.9-GB-SUBS';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/375Dsubs_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% 
% name{2} = 'YSU-3.9-TD-SUBS';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/375Dsubs_edmf0_pbl1_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';

%% BLOCK

vars = {'Times', 'PBLH', 'QVAPOR', 'QCLOUD', 'RTHRATEN', 'T', 'PH', 'PHB', 'EXCH_H', 'U', 'V', 'W', 'P', 'PB', 'DNW', 'LH', 'HFX'}; 


varsMYNN = {'QKE', 'EL_PBL', 'QSHEAR', 'QBUOY', 'QDISS', 'QWT', 'SH3D', 'DTKE'};
% varsMYNN = [];

for n = 1:length(path)
	
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

LESps = '~/Dropbox/UCSD/WRF PBL Parameterization/Code & Data/Data/DYCOMS_ocean_run/cgils_s12_ctl.ps.nc';
LESts = '~/Dropbox/UCSD/WRF PBL Parameterization/Code & Data/Data/DYCOMS_ocean_run/cgils_s12_ctl.ts.nc';

%% Load LES
loadLES;

%% Get height & times
% wrf(1).PBLH = wrf(2).PBLH; % have to put this stupid hack or normalization doesn't work for wrong pblh case
for n = 1:length(wrf)
	wrf(n).T = wrf(n).T + 300;
	wrf(n).time = datenum(wrf(n).Times');
	
	temp=(wrf(n).PH+wrf(n).PHB)./9.81;
	wrf(n).height=(temp(:,:,1:size(temp,3)-1,:)+temp(:,:,2:size(temp,3),:))/2;
	wrf(n).heightstag = temp;
	
	% Normalize height
	for i = 1:size(wrf(1).T, 1)
		for j = 1:size(wrf(1).T, 2)
			
% 			for t = 1:size(wrf(1).T, 4)
% 				wrf(n).heightstag(i,j,1,t) = 0.5 * ( wrf(n).height(i,j,1,t) );
% 				for k = 2:size(wrf(1).T, 3)
% 					wrf(n).heightstag(i,j,k,t) = 0.5 * ( wrf(n).height(i,j,k,t) + wrf(n).height(i,j,k-1,t) );
% 				end
% 				wrf(n).heightstag(i,j,size(wrf(1).T, 3) + 1,t) = 99999;
% 			end
			
% 			for k = 1:size(wrf(n).T, 3)
% 				for t = 2:size(wrf(n).T, 4)
% 					wrf(n).znorm(i,j,k,t) = wrf(n).height(i,j,k,t) ./ wrf(n).PBLH(i, j, t);
% 					wrf(n).znormstag(i,j,k,t) = wrf(n).heightstag(i,j,k,t) ./ wrf(n).PBLH(i, j, t);
% 				end
% 			end
		end
	end
end

%% Load LWP
% Domain averaged? Take single point for now
i = 25; j = 25;

%% PBLH
% 
% % changed all to corrected pblh
% 
figure; plot((1:size(wrf(1).PBLH,3))./4,squeeze(wrf(1).PBLH(25,25,:)), 'r')
hold on; plot((1:size(wrf(2).PBLH,3))./4,squeeze(wrf(2).PBLH(25,25,:)), 'b--')
plot(les.timeTS./3600, les.zi1, 'k')
% plot(squeeze(wrf(1).DELTA_YSU(1,1,:)), 'b-*')
% plot(squeeze(wrf(2).DELTA_YSU(1,1,:)), 'r-o')
legend(name{1}, name{2}, 'LES', 'location' ,'best'); legend('boxoff');
set(gca, 'xtick', 0:6:24)
ylabel('PBLH [m]')
xlabel('Time [hr]')
title PBLH

%% LWP
i = 25; j = 25;
for n = 1:length(wrf)
	for t = 1:size(wrf(n).T, 4)
		wrf(n).lwp(i, j, t) = trapz( squeeze(wrf(n).height( i, j, :, t )), squeeze(wrf(n).QCLOUD( i, j, :, t )) ) .* densityAir;
	end
end

figure;
hold on
for n = 1:length(wrf)
	plot((0:size(wrf(n).lwp,3) - 1)./4,squeeze(wrf(n).lwp(i,j,:)), plotList{n})
end

plot(les.timeTS./3600, les.lwp./1000, 'k')
legend([name, {'LES'}], 'location' ,'best'); legend('boxoff');
set(gca, 'xtick', 0:2:24)
ylabel('LWP [kg m-2]')
xlabel('Time [hr]')
title LWP
xlim([0 24])
hgexport(gcf, ['~/Dropbox/tempPlots/' 'LWP_' name{1} '_' name{2} '.png'], png)

%% TV & TL

for n = 1:length(wrf)
i = 25; j = 25;
			for k = 1:size(wrf(n).T, 3)
				for t = 1:size(wrf(n).T, 4)
					
					wrf(n).TV(i, j, k, t) = (wrf(n).T(i, j, k, t) + 0) .* ( 1 + 0.61 .* wrf(n).QVAPOR(i, j, k, t) - wrf(n).QCLOUD(i, j, k, t) );
					wrf(n).TL(i, j, k, t) = (wrf(n).T(i, j, k, t) + 0) - Lv/cpAir .* ( (wrf(n).P(i, j, k, t) + wrf(n).PB(i, j, k, t)) ./ 100000) .^ (R/cpAir) .* (wrf(n).QCLOUD(i, j, k, t));
					
				end
			end

end


%%
n1 = 1; n2 = 2; % Compare between configurations
i = 25; j = 25;

for n = 1:length(wrf)
	wrf(n).QT = wrf(n).QCLOUD + wrf(n).QVAPOR;% + wrf(n).QRAIN;
end

% plotVars = {'TL', 'QCLOUD', 'SH3D', 'EL_PBL'}; heightVars = {'height', 'height', 'height', 'heightstag'}; x = [285, 310; 0, 0.001; 0, 1.2; 0 500];
% plotVars = {'U', 'V', 'EXCH_H', 'QCLOUD'};
% plotVars = {'T', 'QT', 'QCLOUD', 'U'}; heightVars = {'height', 'height', 'height', 'height'};
% % % % % % % % % % % % % % % plotVars = {'EL_PBL', 'QKE', 'SH3D', 'EXCH_H'};
% % % % % % % % % % % % % % % heightVars = {'heightstag', 'height', 'height', 'heightstag'};

% plotVars = {'T', 'QCLOUD', 'QVAPOR', 'EL_PBL'}; heightVars = {'height', 'height', 'height', 'heightstag'}; x = [285, 310; 0, 0.001; 0, 0.1; 0 500];
% plotVars = {'TL', 'QCLOUD', 'RTHRATEN', 'QKE'}; heightVars = {'height', 'height', 'height', 'height'}; x = [285, 310; 0, 0.001; -0.0005, 0.0005; 0 1.5];
plotVars = {'TL', 'QCLOUD', 'EXCH_H', 'QKE'}; heightVars = {'height', 'height', 'heightstag', 'height'}; x = [285, 310; 0, 0.001; 0, 75; 0 1];
% plotVars = {'TL', 'QCLOUD', 'EXCH_H', 'EL_PBL'}; heightVars = {'height', 'height', 'heightstag', 'heightstag'}; x = [285, 310; 0, 0.001; 0, 300; 0 500];

plotVars = {'TL', 'QCLOUD', 'EXCH_H', 'W'}; heightVars = {'height', 'height', 'heightstag', 'heightstag'}; x = [285, 310; 0, 0.001; 0, 75; -0.005 0.005];

yTop = 3000;
% x = [280, 305; 0, 3000; 0, 2; 0 1500];
% x = [-15, 15; -15, 15; 0, 300; 0 .005];
% x = [285, 310; 0, 0.015; 0, .005; -5 8];
% % % % % % % % % % % % % % x = [0, 200; 0, 1; 0, 1.5; 0 100];

figure('position', [50 50 1400 900])
% for t = 1:480:length(wrf(1).time)
for t = 1:4:length(wrf(1).time)
	for v = 1:length(plotVars)
		
		% Top
		subplot(2,length(plotVars),v)
		plot( squeeze(wrf(n1).(plotVars{v})(i, j, :, t)), squeeze(wrf(n1).(heightVars{v})(i, j, :, t)) );
		hold on
		plot( [-1000 1000], [wrf(n1).PBLH(i,j,t), wrf(n1).PBLH(i,j,t)], '--k')
		hold off
		title(plotVars{v})
		ylim([0 yTop])
		xlim(x(v,:))
		grid on
		
		% Bottom
		subplot(2,length(plotVars),length(plotVars) + v)
		plot( squeeze(wrf(n2).(plotVars{v})(i, j, :, t)), squeeze(wrf(n2).(heightVars{v})(i, j, :, t)) );
		hold on
		plot( [-1000 1000], [wrf(n2).PBLH(i,j,t), wrf(n2).PBLH(i,j,t)], '--k')
		hold off
		title(plotVars{v})
		ylim([0 yTop])
		xlim(x(v,:))
		grid on
		
	end
	text(0.4 .* x(4,2), 0.95.*yTop, ['t = ' num2str(t)])
	text(0.2 .* x(4,2), 0.8.*yTop, name{n1})
	text(0.2 .* x(4,2), 0.6.*yTop, name{n2})
	text(0.2 .* x(4,2), 0.4.*yTop, ['lwp1=' num2str(wrf(n1).lwp(i,j,t))])
	text(0.2 .* x(4,2), 0.3.*yTop, ['lwp2=' num2str(wrf(n2).lwp(i,j,t))])
	hgexport(gcf, ['~/Dropbox/tempPlots/' 'THERMOprofiles_' sprintf('%.2d', t) name{1} '_' name{2} '.png'], png)
	pause()
end

%% LWP Map

wrf(1).XLONG = wrfLoad(path{1}, {'XLONG'});
wrf(1).XLAT = wrfLoad(path{1}, {'XLAT'});

figure('position', [50 50 1400 900])
for t = 1:4:size(wrf(n).T, 4)
	for n = 1:length(wrf)
		for i = 1:size(wrf(1).QCLOUD, 1)
			for j = 1:size(wrf(1).QCLOUD, 2)
				lwp(i, j) = trapz( squeeze(wrf(n).height( i, j, :, t )), squeeze(wrf(n).QCLOUD( i, j, :, t )) ) .* densityAir;
			end
		end
		subplot(1,2,n)
		
		h = pcolor( squeeze(wrf(1).XLONG(:,:,1)), squeeze(wrf(1).XLAT(:,:,1)), lwp);
		hold on
		interval = 4:4:50;
		quiver( squeeze(wrf(1).XLONG(interval, interval,1)), squeeze(wrf(1).XLAT(interval,interval,1)), squeeze(wrf(n).U(interval,interval, 10, t)), squeeze(wrf(n).V(interval,interval, 10, t)), 'r');
		hold off
		set(h, 'edgecolor', 'none')
% 		USBoundaries
		caxis([0 0.3])
		xlabel('Longitude [deg]'), ylabel('Latitude [deg]'), title('LWP')
	end
	title(datestr(wrf(1).time(t)))
	colorbar
	hgexport(gcf, ['~/Dropbox/tempPlots/' 'mapLWP_' sprintf('%.3d', t) name{1} '_' name{2} '.png'], png)
	pause(0.1)
	clf
end

return

% hgexport(gcf, ['tempPlots/' 'LWP_' name{1} '_' name{2} '.png'], png)

%% TKE BUDGET
n1 = 1; n2 = 2; % Compare between configurations
i = 25; j = 25;

yTop = 2000;

figure('position', [50 50 1400 900])
% for t = 1:480:length(wrf(1).time)
for t = 1:1:length(wrf(1).time)
	for v = 1:length(plotVars)
		% Top
		clf
		hold on
		plot(squeeze(wrf(n1).QSHEAR(i, j, :, t)), squeeze(wrf(n1).heightstag(i, j, :, t)), '--g')
		plot(squeeze(wrf(n1).QBUOY(i, j, :, t)), squeeze(wrf(n1).heightstag(i, j, :, t)), 'y--')
		plot(squeeze(wrf(n1).QDISS(i, j, :, t)), squeeze(wrf(n1).heightstag(i, j, :, t)), 'r--')
		plot(squeeze(wrf(n1).QWT(i, j, :, t)), squeeze(wrf(n1).heightstag(i, j, :, t)), 'b--')
		plot(squeeze(wrf(n1).DTKE(i, j, :, t)), squeeze(wrf(n1).height(i, j, :, t)), 'k--')
		hold off
		
		ylim([0 yTop])
		title(num2str(t))
		grid on
		
	end
	% 	hgexport(gcf, ['~/Dropbox/Presentations/SCM-VertAdv/' '_profile_' sprintf('%.2d', t) name1 '_' name2 '.png'], png)
	pause()
end

%% Four panel TKE budget
n1 = 1; n2 = 2;
hoursEnd = 1:4;
yTop = 1200;
x1 = -0.03; x2 = 0.03;

figure('position', [50 50 1400 900])
for h = 1:length(hoursEnd)
	
	resolution = 15; % Minutes
	endIdx = 1 + 60 ./ resolution .* hoursEnd(h); startIdx = endIdx - 60 ./ resolution; 
	
	subplot(2,length(hoursEnd),h)
	plot( mean( squeeze( wrf(n1).QBUOY(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n1).heightstag(i, j, :, startIdx:endIdx) ), 2), 'color', [1 0.6 0])
	hold on
	plot( mean( squeeze( wrf(n1).QSHEAR(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n1).heightstag(i, j, :, startIdx:endIdx) ), 2), 'g--')
	plot( mean( squeeze( wrf(n1).QDISS(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n1).heightstag(i, j, :, startIdx:endIdx) ), 2), 'r--')
	plot( mean( squeeze( wrf(n1).QWT(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n1).heightstag(i, j, :, startIdx:endIdx) ), 2), 'b--')
	
	ylim([0 yTop])
	xlim([x1 x2])
	xlabel('m2 s-2')
	
	subplot(2,length(hoursEnd),length(hoursEnd) + h)
	plot( mean( squeeze( wrf(n2).QBUOY(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n2).heightstag(i, j, :, startIdx:endIdx) ), 2), 'color', [1 0.6 0])
	hold on
	plot( mean( squeeze( wrf(n2).QSHEAR(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n2).heightstag(i, j, :, startIdx:endIdx) ), 2), 'g--')
	plot( mean( squeeze( wrf(n2).QDISS(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n2).heightstag(i, j, :, startIdx:endIdx) ), 2), 'r--')
	plot( mean( squeeze( wrf(n2).QWT(i, j, :, startIdx:endIdx) ), 2), mean( squeeze( wrf(n2).heightstag(i, j, :, startIdx:endIdx) ), 2), 'b--')
	title(['t = ' num2str(h) ' hr'])
	ylim([0 yTop])
	xlim([x1 x2])
	
	hgexport(gcf, ['~/Dropbox/tempPlots/' 'TKE_' sprintf('%.2d', t) name{1} '_' name{2} '.png'], png)
	
end

%% Compute BL means
for n = length(wrf):-1:1
	for i = 25% 1:size(wrf(1).T, 1)
		for j = 25% 1:size(wrf(1).T, 2)
			for t = 1:size(wrf(1).T, 4)
				wrf(n).kpbl(i, j, t) = findMatch( squeeze(wrf(n).height( i, j, :, t)), wrf(n).PBLH(i, j, t) ) - 1;
				% determine cloud base
				wrf(n).zb(i,j,t) =  wrf(n).height(i, j, find(squeeze(wrf(n).QCLOUD(1,1,:,10)),  1, 'first'), t);
% 				if n < 3
% 					wrf(1).kpbl(i, j, t) = wrf(2).kpbl(i, j, t); % Force original ACM2 to have a correct PBL height
% % 					wrf(2).kpbl(i, j, t) = wrf(3).kpbl(i, j, t); % Force original ACM2 to have a correct PBL height
% 				end
% 				if n <= 3
% 					wrf(1).kpbl(i,j,t) = 32;
% 					wrf(2).kpbl(i,j,t) = 32;
% 				end
				wrf(n).meanThetaL(i, j, t) = sum( squeeze(wrf(n).TL( i, j, 1:wrf(n).kpbl( i, j, t ), t)) .* squeeze( wrf(n).DNW( 1:wrf(n).kpbl( i, j, t ), t ) ) ) ./ sum( wrf(n).DNW( 1:wrf(n).kpbl( i, j, t ), t ) );
				wrf(n).meanQT(i, j, t) = sum( squeeze(wrf(n).QVAPOR(i, j, 1:wrf(n).kpbl( i, j, t), t) + wrf(n).QCLOUD( i, j, 1:wrf(n).kpbl( i, j, t ), t)) .* squeeze( wrf(n).DNW( 1:wrf(n).kpbl( i, j, t ), t ) ) ) ./ sum( wrf(n).DNW( 1:wrf(n).kpbl( i, j, t ), t ) );
			end
		end
	end
end

%% Normalize and time average

zRange = 0:0.05:2;

n = 1; % n of MYNN case
i = 25; j = 25;

zinv = 800;

for t = 2:size(wrf(1).T, 4)
	haveDTKE = squeeze(wrf(n).DTKE(i, j, :, t));
	haveQKE = squeeze(wrf(n).QKE(i, j, :, t));
	haveEL = squeeze(wrf(n).EL_PBL(i, j, :, t));
	haveS = squeeze(wrf(n).QSHEAR(i, j, :, t)) ./ haveEL;
	haveB = squeeze(wrf(n).QBUOY(i, j, :, t)) ./ haveEL;
	haveD = squeeze(wrf(n).QDISS(i, j, :, t)) ./ haveEL;
	haveT = squeeze(wrf(n).QWT(i, j, :, t)) ./ haveEL;
	
% 	haveTSQ = squeeze(wrf(n).TSQ(i, j, :, t));
% 	haveQSQ = squeeze(wrf(n).QSQ(i, j, :, t));
% 	haveCOV = squeeze(wrf(n).COV(i, j, :, t));
% 	haveWSQ = squeeze(wrf(n).WSQ(i, j, :, t));
	
	haveSH3D = squeeze(wrf(n).SH3D(i, j, :, t));
	haveEXCH_H = squeeze(wrf(n).EXCH_H(i, j, :, t));
	
% 	haveZ = squeeze(wrf(n).height(i,j,:,t) ./ wrf(n).PBLH(i, j, t));
% 	haveZstag = squeeze(wrf(n).heightstag(i,j,:,t) ./ wrf(n).PBLH(i, j, t));

	haveZ = squeeze(wrf(n).height(i,j,:,t) ./ zinv);
	haveZstag = squeeze(wrf(n).heightstag(i,j,:,t) ./ zinv);
	
	wrf(n).iDTKE(i, j, :, t) = interp1(haveZ, haveDTKE, zRange);
	wrf(n).iQKE(i, j, :, t) = interp1(haveZ, haveQKE, zRange);
	wrf(n).iS(i, j, :, t) = interp1(haveZstag, haveS, zRange);
	wrf(n).iB(i, j, :, t) = interp1(haveZstag, haveB, zRange);
	wrf(n).iD(i, j, :, t) = interp1(haveZstag, haveD, zRange);
	wrf(n).iT(i, j, :, t) = interp1(haveZstag, haveT, zRange);
	wrf(n).iEL(i, j, :, t) = interp1(haveZstag, haveEL, zRange);
	
% 	wrf(n).iTSQ(i, j, :, t) = interp1(haveZ, haveTSQ, zRange);
% 	wrf(n).iQSQ(i, j, :, t) = interp1(haveZ, haveQSQ, zRange);
% 	wrf(n).iCOV(i, j, :, t) = interp1(haveZ, haveCOV, zRange);
% 	wrf(n).iWSQ(i, j, :, t) = interp1(haveZ, haveWSQ, zRange);
	
	wrf(n).iSH3D(i, j, :, t) = interp1(haveZ, haveSH3D, zRange);
	wrf(n).iEXCH_H(i, j, :, t) = interp1(haveZstag, haveEXCH_H, zRange);
	
end


% 2 - 4
tRange =  (1+(2*4)):(4*4 + 1); % 4-6
lesRange = 2:4;
% 
% 11 - 13
% tRange =  (1+(11*60)):(13*60 + 1); % 11-13
% lesRange = 11:13;


for z = 1:length(zRange)
	avgDTKE(n,z) = mean(wrf(n).iDTKE(i,j,z,tRange), 4);
	avgQKE(n,z) = mean(wrf(n).iQKE(i,j,z,tRange), 4);
	avgS(n,z) = mean(wrf(n).iS(i,j,z,tRange), 4);
	avgB(n,z) = mean(wrf(n).iB(i,j,z,tRange), 4);
	avgD(n,z) = mean(wrf(n).iD(i,j,z,tRange), 4);
	avgT(n,z) = mean(wrf(n).iT(i,j,z,tRange), 4);
	avgEL(n,z) = mean(wrf(n).iEL(i,j,z,tRange), 4);
	
% 	avgTSQ(n,z) = mean(wrf(n).iTSQ(i,j,z,tRange), 4);
% 	avgQSQ(n,z) = mean(wrf(n).iQSQ(i,j,z,tRange), 4);
% 	avgCOV(n,z) = mean(wrf(n).iCOV(i,j,z,tRange), 4);
% 	avgWSQ(n,z) = mean(wrf(n).iWSQ(i,j,z,tRange), 4);
% 	avgWSQms(n,z) = mean(wrf(n).iWSQ(i,j,z,tRange) .* wrf(n).iQKE(i,j,z,tRange), 4);
	
	avgSH3D(n,z) = mean(wrf(n).iSH3D(i,j,z,tRange), 4);
	avgEXCH_H(n,z) = mean(wrf(n).iEXCH_H(i,j,z,tRange), 4);
	
end

avgZB(n) = mean(wrf(n).zb(i,j,tRange));

% Average LES
startIdx = 1 + lesRange(1)*3600/20/60; endIdx = 1 + lesRange(end)*3600/20/60;
startIdx2 = 1 + lesRange(1)*3600/20; endIdx2 = 1 + lesRange(end)*3600/20;
midIdx = floor( (startIdx2 + endIdx2) / 2);
for z = 1:size(les.wthl,1)
	les.avgTKE(z) = mean(les.tke(z,startIdx:endIdx));
	les.avgWTHL(z) = mean(les.wthl(z,startIdx:endIdx));
	les.avgWQT(z) = mean(les.wqt(z,startIdx:endIdx));
	les.avgWW(z) = mean(les.ww(z,startIdx:endIdx));
	les.avgBUOY(z) = mean(les.buoy(z,startIdx:endIdx));
	les.avgSHEAR(z) = mean(les.shr(z,startIdx:endIdx));
	les.avgTRANS(z) = mean(les.trans(z,startIdx:endIdx));
	les.avgDISS(z) = mean(les.diss(z,startIdx:endIdx));
end
les.avgZC = mean(les.zc(startIdx2:endIdx2));
les.avgZB = mean(les.zb(startIdx2:endIdx2));
norm = 1 ./ les.zi1(floor( (startIdx2 + endIdx2) / 2 ));


%% Plot
%%%%%%%%%%% 2-hr avg %%%%%%%%%%%
figure('position', [50 50 1200 600])

subplot(1,3,1)
plot( .5 .* avgQKE(n,:), zRange, 'k' ) % 2*TKE
hold on
plot(les.avgTKE(2:end), les.zt(2:end) .* norm, 'k--')

plot([-1000 1000], [1, 1], 'k:', 'linewidth', 1) % PBLH
plot([-1000 1000], [les.avgZB.*norm, les.avgZB.*norm], 'k:', 'linewidth', 1) % PBLH
plot([-1000 1000], [avgZB(1)./zinv, avgZB(1)./zinv], 'r-.', 'linewidth', 1) % zb MYNN
plot([0 0], [0 10000], 'k:') % 0 origin line
legend('TKE', 'TKE_{les}')

xlim([0 1])
ylim([0 1.5])
xlabel('TKE [m2 s-2]')
ylabel('z / z_i [-]')
title TKE

subplot(1,3,2)
% plot( 100 .* avgDTKE(n,:), zRange, 'k' ) % DTKE
hold on
% plot( 100.* squeeze(wrf(n).TKE_PBL(i, j, :, t)), squeeze(wrf(n).heightstag(i, j, :, t)), 'k^--' ) % TKE FROM PBL -- VERIFIED MATCH
plot( 100.* avgS(n,:), zRange, 'g' ) % Shear production
plot( 100.* avgB(n,:), zRange, 'm' ) % Buoyant production
plot(100.*les.avgSHEAR(1:2:end), les.zt(1:2:end) .* norm, 'go', 'markersize', 2)
plot(100.*les.avgBUOY(2:2:end), les.zm(2:2:end) .* norm, 'mo', 'markersize', 2)

% divide by L

plot([-1000 1000], [1, 1], 'k:', 'linewidth', 1) % PBLH
% plot([-1000 1000], [les.avgZC.*norm, les.avgZC.*norm], 'k', 'linewidth', 1) % PBLH
plot([-1000 1000], [les.avgZB.*norm, les.avgZB.*norm], 'k:', 'linewidth', 1) % zb LES
plot([-1000 1000], [avgZB(1)./zinv, avgZB(1)./zinv], 'r-.', 'linewidth', 1) % zb MYNN
plot([0 0], [0 10000], 'k:') % 0 origin line
legend('S', 'B', 'S_{les}', 'B_{les}')

ylim([0 1.5])
xlim([-0.1 0.45])
xlabel('[cm2 s-3]')
title Production

subplot(1,3,3)

plot(-100.* avgD(n,:), zRange, 'r' ) % Dissipation
hold on
plot(100.* avgT(n,:), zRange, 'b' ) % Transport
plot(-100.*les.avgDISS(1:2:end), les.zm(1:2:end) .* norm, 'ro', 'markersize', 2)
plot(100.*les.avgTRANS(1:2:end), les.zm(1:2:end) .* norm, 'bo', 'markersize', 2)

plot([-1000 1000], [1, 1], 'k:', 'linewidth', 1) % PBLH
plot([-1000 1000], [les.avgZB.*norm, les.avgZB.*norm], 'k:', 'linewidth', 1) % PBLH
plot([-1000 1000], [avgZB(1)./zinv, avgZB(1)./zinv], 'r-.', 'linewidth', 1) % zb MYNN
plot([0 0], [0 10000], 'k:') % 0 origin line
legend('D', 'T', 'D_{les}', 'T_{les}', 'location', 'southwest')

ylim([0 1.5])
xlim([-0.7 0.4])
xlabel('[cm2 s-3]')
title ('Dissip. & Transport')

%% INT TKE

n = 1;
i = 25; j = 25;
for n = 1:length(wrf)
	for t = 1:length(wrf(n).time)
		% 	wrf(n).intTKE(1, 1, t) = sum( squeeze(wrf(n).TKE_PBL( 1, 1, 1:70, t)) .* squeeze( wrf(n).DNW( 1:70, t ) ) ) ./ sum( wrf(n).DNW( 1:70, t ) );
		wrf(n).intTKE(i, j, t) = trapz( squeeze(wrf(n).height( i, j, :, t )), 0.5 .* squeeze(wrf(n).QKE( i, j, :, t )) );
	end
end
% for t = 1:length(les.time)
% 	les.TKEintManual(t) = trapz(les.zt,les.tke(:,t));
% end

figure; 
plot( (1:length(wrf(1).time)) ./ 4, squeeze(wrf(1).intTKE(i,j,:)), 'r'); hold on
plot( (1:length(wrf(2).time)) ./ 4, squeeze(wrf(2).intTKE(i,j,:)), 'b');
% plot( (1:length(wrf(3).time)) ./ 60, squeeze(wrf(3).intTKE(i,j,:)), 'g');
hold on
plot(les.timeTS./3600, les.tkeint, 'k')
% plot(les.time./3600, les.TKEintManual, 'r--')
% legend('MYNN', 'LES', 'location', 'best')
h = legend(name{1}, name{2}, 'LES', 'location', 'best');
legend('boxoff');
set(gca, 'xtick', 0:2:24)
xlim([0 24])
grid on
ylabel('Vert. int. TKE [m3 s-2]')
xlabel('Time [hr]')

hgexport(gcf, ['~/Dropbox/tempPlots/' 'INTTKE_' name{1} '_' name{2} '.png'], png)

%% Algorithm to determine PBLH based on Cloud Top height

for n = 1:length(wrf)
	i = 25; j = 25;
	for t = 1:size(wrf(n).T, 4)
		if ~isempty(find(squeeze(wrf(n).QCLOUD(i, j, :, t)), 1, 'last'))
			wrf(n).ctop(i, j, t) = find(squeeze(wrf(n).QCLOUD(i, j, :, t)), 1, 'last');
		else
			wrf(n).ctop(i, j, t) = 0;
		end
	end
end


%% Compute BL means
for n = 1:length(wrf)
	wrf(n).DNW = ncread(path{n}, 'DNW');
end

for n = length(wrf):-1:1
	i = 25; j = 25;
	for t = 1:size(wrf(1).T, 4)
		if wrf(n).ctop(i, j, t) ~= 0
			% 				wrf(n).kpbl(i, j, t) = findMatch( squeeze(wrf(n).height( i, j, :, t)), wrf(n).PBLH(i, j, t) ) - 1;
% 			wrf(n).kpbl(i, j, t) = findMatch( squeeze(wrf(n).height( i, j, :, t)), wrf(n).ctop(i,j,t) ) - 1;
			wrf(n).kpbl(i, j, t) = wrf(n).ctop(i,j,t) - 1;
			% determine cloud base
			wrf(n).zb(i,j,t) =  wrf(n).height(i, j, find(squeeze(wrf(2).QCLOUD(1,1,:,10)),  1, 'first'), t);
			% 				if n < 3
			% 					wrf(1).kpbl(i, j, t) = wrf(2).kpbl(i, j, t); % Force original ACM2 to have a correct PBL height
			% % 					wrf(2).kpbl(i, j, t) = wrf(3).kpbl(i, j, t); % Force original ACM2 to have a correct PBL height
			% 				end
			% 				if n <= 3
			% 					wrf(1).kpbl(i,j,t) = 32;
			% 					wrf(2).kpbl(i,j,t) = 32;
			% 				end
			
			wrf(n).meanThetaL(i, j, t) = sum( squeeze(wrf(n).TL( i, j, 1:wrf(n).kpbl( i, j, t ), t)) .* squeeze( wrf(n).DNW( 1:wrf(n).kpbl( i, j, t ), t ) ) ) ./ sum( wrf(n).DNW( 1:wrf(n).kpbl( i, j, t ), t ) );
			wrf(n).meanQT(i, j, t) = sum( squeeze(wrf(n).QVAPOR(i, j, 1:wrf(n).kpbl( i, j, t), t) + wrf(n).QCLOUD( i, j, 1:wrf(n).kpbl( i, j, t ), t)) .* squeeze( wrf(n).DNW( 1:wrf(n).kpbl( i, j, t ), t ) ) ) ./ sum( wrf(n).DNW( 1:wrf(n).kpbl( i, j, t ), t ) );
		else
			wrf(n).kpbl(i,j,t) = 1;
			wrf(n).zb(i,j,t) = 1;
			wrf(n).meanThetaL(i, j, t) = 0;
			wrf(n).meanQT(i, j, t) = 0;
		end
	end
end


%% BL Means

% THETA_L
figure('position', [50 50 600 300]);
plot((1:length(wrf(1).time))./4,squeeze(wrf(1).meanThetaL(i,j,:)), 'r')
hold on; plot((1:length(wrf(2).time))./4,squeeze(wrf(2).meanThetaL(i,j,:)), 'b--')
% plot((1:length(wrf(3).time))./4,squeeze(wrf(3).meanThetaL(i,j,:)), 'g')
plot(les.time./3600, les.meanTL, 'k')
% plot(squeeze(wrf(1).DELTA_YSU(1,1,:)), 'b-*')
% plot(squeeze(wrf(2).DELTA_YSU(1,1,:)), 'r-o')
legend(name{1}, name{2}, 'LES', 'location' ,'best'); legend('boxoff');
set(gca, 'xtick', [0:6:24])
ylabel('\theta_l [K]')
xlabel('Time [hr]')
title <\theta_l>_{BL}

% hgexport(gcf, [plotDir type '_meanTL_' name1 '_' name2 '_' name3 '.png'], png)

% QT
figure('position', [50 50 600 300]);
plot((1:length(wrf(1).time))./4,squeeze(wrf(1).meanQT(i,j,:)), 'r')
hold on; plot((1:length(wrf(2).time))./4,squeeze(wrf(2).meanQT(i,j,:)), 'b--')
% plot((1:length(wrf(1).time))./4,squeeze(wrf(3).meanQT(i,j,:)), 'g')
plot(les.time./3600, les.meanQT.*1e-3, 'k')
% plot(squeeze(wrf(1).DELTA_YSU(1,1,:)), 'b-*')
% plot(squeeze(wrf(2).DELTA_YSU(1,1,:)), 'r-o')
legend(name{1}, name{2}, 'LES', 'location' ,'best'); legend('boxoff');
set(gca, 'xtick', [0:6:24])
ylabel('q_T [kg kg-1]')
xlabel('Time [hr]')
title <q_T>_{BL}

% hgexport(gcf, [plotDir type '_meanQT_' name1 '_' name2 '_' name3 '.png'], png)
