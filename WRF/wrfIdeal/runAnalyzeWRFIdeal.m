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
name{11} = 'MYNN-3.9-EDMF';
path{11} = [dycomsrf01 'edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% 
% % 0.5 EL
name{6} = 'MYNN-3.9-0.5EL';
path{6} = [dycomsrf01 '0.5EL_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% 0.5 EL PROFILE
name{2} = 'MYNN-3.9-0.5ELPROF';
path{2} = [dycomsrf01 '0.5ELPROFILE_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % 
% % 
% 2.0 EL
name{5} = 'MYNN-3.9-1.5EL';
path{5} = [dycomsrf01 '1.5_EL_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% 2.0 EL PROFILE
name{3} = 'MYNN-3.9-1.5ELPROF';
path{3} = [dycomsrf01 '1.5_ELPROFILE_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % 
% % % 0.5 EL
% % 
name{7} = 'MYNN-3.9-EDMF-0.5EL';
path{7} = [dycomsrf01 '0.5EL_edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% 0.5 EL PROFILE
name{8} = 'MYNN-3.9-EDMF-0.5ELPROF';
path{8} = [dycomsrf01 '0.5ELPROFILE_edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % 
% % 
% % 2.0 EL
name{9} = 'MYNN-3.9-EDMF-1.5EL';
path{9} = [dycomsrf01 '1.5_EL_edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% 2.0 EL PROFILE
name{10} = 'MYNN-3.9-EDMF-1.5ELPROF';
path{10} = [dycomsrf01 '1.5_ELPROFILE_edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];

% YSU and RAD5 MYNN experiment

name{4} = 'YSU-3.9-TD';
path{4} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/375Dsubs_edmf0_pbl1_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';

% name{12} = 'MYNN-3.9-GB-RAD5';
% path{12} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/edmf0_pbl5_sw5_lw5_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';

% name{13} = 'DK';
% path{13}= '/mnt/lab_48tb2b/users/dksahu/WRF-DA/20130603/WRFDAintmt-tune1/wrfout_d02_2013-06-03_12_00_00';

%% BUOY 3D sims

clear path
path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfSCM/scmPaperFINAL/3D/ysutd/wrfout_d01_2015-07-10_08_00_00';
name{1} = 'YSU-TD';

path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfSCM/scmPaperFINAL/3D/ysubuoy/wrfout_d01_2015-07-10_08_00_00';
name{2} = 'YSU-BUOY';

path{3} = '/mnt/lab_48tb1/users/hyang/database/wrfSCM/scmPaperFINAL/3D/acm2/wrfout_d01_2015-07-10_08_00_00';
name{3} = 'ACM2';

path{4} = '/mnt/lab_48tb1/users/hyang/database/wrfSCM/scmPaperFINAL/3D/acm2buoy/wrfout_d01_2015-07-10_08_00_00';
name{4} = 'ACM2-BUOY';

path{5} = '/mnt/lab_48tb1/users/hyang/database/wrfSCM/scmPaperFINAL/3D/mynn/wrfout_d01_2015-07-10_08_00_00';
name{5} = 'MYNN';

path{6} = '/mnt/lab_48tb1/users/hyang/database/wrfSCM/scmPaperFINAL/3D/acm2buoyA15/wrfout_d01_2015-07-10_08_00_00';
name{6} = 'ACM2-BUOY-A15';
% 
% path{7} = '/mnt/lab_48tb1/users/hyang/database/wrfSCM/scmPaperFINAL/3D/acm2buoyA30/wrfout_d01_2015-07-10_08_00_00';
% name{7} = 'ACM2-BUOY-A30';

for n = 1:length(path)
	try
	out = analyzeWRFIdeal(path{n}, name{n});
	catch
	end
end

%%
% 
% name{1} = 'MYNN-3.6';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.6/wrfout_d01_20150710_0800.nc';
% % % name{2} = 'MYNN-3.6-NBL';
% % % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.6noBouLac/wrfout_d01_20150710_0800.nc';
% % % name{2} = 'MYNN-3.6-MYNN3';
% % % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.6MYNN3/wrfout_d01_20150710_0800.nc';
% % % % % % name{1} = 'MYNN-3.7.1';
% % % % % % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1mynn/wrfout_d01_20150710_0800.nc';
% % % % % % name{2} = 'MYNN-3.7.1-0F';
% % % % % % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1_equator/wrfout_d01_20150710_0800.nc';
% 
% % name{1} = 'MYNN-3.7.1';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1mynn_0f/wrfout_d01_20150710_0800.nc';
% % name{1} = 'MYNN-3.8.1-ML';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1mynnML2_0f/wrfout_d01_20150710_0800.nc';
% % name{2} = 'MYNN-3.8.1-0f';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1mynn_0f/wrfout_d01_20150710_0800.nc';
% % 
% % 
% % name{1} = 'YSU-TD-3.7.1';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1ysutd_0f/wrfout_d01_20150710_0800.nc';
% % name{2} = 'YSU-BUOY-3.7.1';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1ysubuoy_0f/wrfout_d01_20150710_0800.nc';
% % name{3} = 'YSU-TD-3.8.1';
% % path{3} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1ysutd_0f/wrfout_d01_20150710_0800.nc';
% 
% 
% % % % name{2} = 'MYNN-3.8.1';
% % % % % % % % % % % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1NEWWIND/wrfout_d01_20150710_0800.nc';
% % % % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1/wrfout_d01_20150710_0800.nc';
% % name{1} = 'MYNN-3.8.1';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1_24hr/wrfout_d01_20150710_0800.nc';
% % name{1} = 'MYNN-3.8.1-TEMF';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1temf/wrfout_d01_20150710_0800.nc';
% 
% % name{2} = 'YSU-3.6';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.6ysu/wrfout_d01_20150710_0800.nc';
% % name{1} = 'YSU-TD-3.8.1';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1ysuTD/wrfout_d01_20150710_0800.nc';
% % name{1} = 'YSU-TD-3.7.1';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1ysuTD/wrfout_d01_20150710_0800.nc';
% % % name{2} = 'YSU-BY-3.7.1';
% % % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1ysubuoy/wrfout_d01_20150710_0800.nc';
% % name{2} = 'YSU-3.8.1';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1ysu/wrfout_d01_20150710_0800.nc';
% 
% %% HY test
% % name{1} = 'MYNN-3.9-GB';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01GB/wrfout_d01_20150710_0800.nc';
% % 
% name{1} = 'MYNN-3.9-GB-51';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfSCM/testGB_39vanillathl_51/wrfout_d01_2015-07-10_08:00:00';
% 
% % 
% % % name{2} = 'MYNN-3.8.1';
% % % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1mynn_0f/wrfout_d01_20150710_0800.nc';
% % name{2} = 'MYNN-3.8.1-51';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.8.1mynn_51/wrfout_d01_20150710_0800.nc';
% % 
% % name{2} = 'MYNN-3.9';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn/wrfout_d01_20150710_0800.nc';
% 
% 
% %%%%%%%%%%%%%%%
% 
% % % Rad 0
% % name{1} = 'MYNN-3.9-GB-0r';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_0rad/wrfout_d01_2015-07-10_08_00_00';
% % name{2} = 'MYNN-3.9-0r';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51_norad/wrfout_d01_20150710_0800.nc';
% % 
% % % Rad 5
% name{1} = 'MYNN-3.9-GB-5r';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_5rad/wrfout_d01_2015-07-10_08_00_00';
% name{2} = 'MYNN-3.9-GB-4r';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_4rad/wrfout_d01_2015-07-10_08_00_00';
% 
% %%
% % % Rad 5
% % name{1} = 'MYNN-3.9-GB-4r-1icloud';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_4rad/wrfout_d01_2015-07-10_08_00_00';
% name{1} = 'MYNN-3.9-GB-4r-0icloud';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_4rad_0icloudbl/wrfout_d01_2015-07-10_08_00_00';
% 
% % ELB
% % name{1} = 'MYNN-3.9-GB-200ELP';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/200ELPROFILE_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% % name{1} = 'MYNN-3.9-GB-2ELB';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/2ELB_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% % 
% % name{1} = 'MYNN-3.9-GB-EDMF1';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% name{2} = 'MYNN-3.9-GB-10ELPROFILE';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/100ELPROFILE_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00'; % actually 50 EL profile
% 
% % 2x mixing (EL and K)
% % name{1} = 'MYNN-3.9-GB-2EL';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/2EL_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% % name{2} = 'MYNN-3.9-GB-0.5EL';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/0.5EL_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% % name{1} = 'MYNN-3.9-GB-4r-0icloud';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_4rad_0icloudbl/wrfout_d01_2015-07-10_08_00_00';
% % name{2} = 'MYNN-3.9-GB-2K';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/2K_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% 
% 
% % name{2} = 'MYNN-3.9-GB-4r-1icloud-1ml';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_4rad_1icloudbl_1ml/wrfout_d01_2015-07-10_08_00_00';
% % name{2} = 'MYNN-3.9-GB-5r-1icloud-1ml';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_5rad_1icloudbl_1ml/wrfout_d01_2015-07-10_08_00_00';
% 
% % name{2} = 'MYNN-3.9-5r';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51_5rad/wrfout_d01_2015-07-10_08_00_00';
% % 
% % % 3.9 no GB
% % name{1} = 'MYNN-3.9-0r';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51_norad/wrfout_d01_20150710_0800.nc';
% % name{2} = 'MYNN-3.9-5r';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51_5rad/wrfout_d01_2015-07-10_08_00_00';
% % 
% % name{1} = 'MYNN-3.9-GB-0r';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_0rad/wrfout_d01_2015-07-10_08_00_00';
% % name{2} = 'MYNN-3.9-GB-5r';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_5rad/wrfout_d01_2015-07-10_08_00_00';
% 
% %% BLOCK DATA SET DYCOMS RF01 EXPERIMENT
% 
% dycomsrf01 = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/';
% 
% % % Control 3.6
% % name{1} = 'MYNN-3.6';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.6/wrfout_d01_20150710_0800.nc';
% % 
% % % Control 3.7.1
% % name{1} = 'MYNN-3.7.1';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.7.1mynn_0f/wrfout_d01_20150710_0800.nc';
% 
% % CONTROL CASE NO EDMF
% name{1} = 'MYNN-3.9';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/3.9mynn_51gb_4rad_0icloudbl/wrfout_d01_2015-07-10_08_00_00'; % CONTROL
% % 
% % % CONTROL CASE WITH EDMF = 1
% % name{2} = 'MYNN-3.9-EDMF';
% % path{2} = [dycomsrf01 'edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % 
% % % 0.5 EL
% % name{1} = 'MYNN-3.9-0.5EL';
% % path{1} = [dycomsrf01 '0.5EL_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % % 0.5 EL PROFILE
% % name{2} = 'MYNN-3.9-0.5ELPROF';
% % path{2} = [dycomsrf01 '0.5ELPROFILE_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % % 
% % % 
% % % 2.0 EL
% % name{1} = 'MYNN-3.9-1.5EL';
% % path{1} = [dycomsrf01 '1.5_EL_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % % 2.0 EL PROFILE
% % name{2} = 'MYNN-3.9-1.5ELPROF';
% % path{2} = [dycomsrf01 '1.5_ELPROFILE_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % % 
% % % % 0.5 EL
% % % 
% % name{1} = 'MYNN-3.9-EDMF-0.5EL';
% % path{1} = [dycomsrf01 '0.5EL_edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % % 0.5 EL PROFILE
% % name{2} = 'MYNN-3.9-EDMF-0.5ELPROF';
% % path{2} = [dycomsrf01 '0.5ELPROFILE_edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % % 
% % % 
% % % 2.0 EL
% % name{1} = 'MYNN-3.9-EDMF-1.5EL';
% % path{1} = [dycomsrf01 '1.5_EL_edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% % % 2.0 EL PROFILE
% % name{2} = 'MYNN-3.9-EDMF-1.5ELPROF';
% % path{2} = [dycomsrf01 '1.5_ELPROFILE_edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2' '/wrfout_d01_2015-07-10_08_00_00'];
% 
% % NOBDY
% % name{1} = 'YSU-3.9-TD-0BDY';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/nobdy_edmf0_pbl1_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% 
% % name{1} = 'MYNN-3.9-GB-0BDY';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/nobdy_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% 
% % WITH SUBS
% 
% % name{1} = 'YSU-3.9-TD-0BDY';
% % path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/375Dsubs_edmf0_pbl1_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% 
% % name{2} = 'MYNN-3.9-GB-SUBS';
% % path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/375Dsubs_edmf0_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% 
% name{2} = 'YSU-3.9-TD-SUBS';
% path{2} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/375Dsubs_edmf0_pbl1_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
% 
% name{1} = 'MYNN-3.9-GB-RAD5';
% path{1} = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/dycomsrf01/results/edmf0_pbl5_sw5_lw5_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00';
