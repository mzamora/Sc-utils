%% Generates "ideal" WRF domain within em_real

% Constants
defineConstants
% Path to raw met_em files

path = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/';

% readFile1 = [path '/dycomsrf01/met_em.d01.2015-07-10_08:00:00.nc'];
% writeFile1 = [path '/dycomsrf01/edited/met_em.d01.2015-07-10_08:00:00.nc'];

% CBL
% readFile = [path '/rap13km/raw/met_em.d01.2015-07-10_08:00:00.nc'];
% writeFile = [path '/rap13km/TRASHrf01/met_em.d01.2015-07-10_08:00:00.nc'];
% writeFile = [path '/rap13km/cbl/met_em.d01.2015-07-10_08:00:00.nc'];

% Path to LES files

% RF01
% LESps = '~/Dropbox/UCSD/WRF PBL Parameterization/Code & Data/Data/DYCOMS_ocean_run/cgils_s12_ctl.ps.nc';
% LESts = '~/Dropbox/UCSD/WRF PBL Parameterization/Code & Data/Data/DYCOMS_ocean_run/cgils_s12_ctl.ts.nc';

LESps = '~/Dropbox/lesData/DYCOMS_ocean_run/cgils_s12_ctl.ps.nc';
LESts = '~/Dropbox/lesData/DYCOMS_ocean_run/cgils_s12_ctl.ts.nc';

% CBL
% LESps = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/rap13km/cbl/dcbl_x.ps.nc';

% Should take data from DYCOMS time period due to radiation forcing

% CGILS S12 CTL
readFile = [path '/rap13km/raw/met_em.d01.2015-07-10_08:00:00.nc'];
writeFile = [path '/rap13km/cgils_s12/met_em.d01.2015-07-10_08:00:00.nc'];
LESps = '/mnt/lab_45d1/database/Sc_group/uclales_output/cgils_s12_ctl/hr_1_3_60min/cgils_s12_ctl.ps.nc';
LESts = '/mnt/lab_45d1/database/Sc_group/uclales_output/cgils_s12_ctl/hr_1_3_60min/cgils_s12_ctl.ts.nc';

%% Notes: met_em approach

% 1) Remake netCDF with same structure
% 2) Redefine vertical coordinates, i.e. num_metgrid_levels = 216
% 3) Replace variables with LES
% 	- Times, PRES, GHT, SKINTEMP, PSFC, RH, (QVAPOR), (QCLOUD), VV, UU, TT, PMSL (SAME AS PSFC)
%   - time, p, zt, tsrf or tskinav, RH, (q - l), (l), v, u, t (theta_l, needs adjustment), PSFC/PMSL = same as in sounding input
% 4) Modify all land values
% 	- VAR, OA1, OL1, etc. are used only for gravitational wave option, so are uneccessary--therefore can 0

% Vars with time dimension
% Everything has Time as final dimension

%% Notes: geo_em approach

% 1) One copy of normal geo
% 2) One copy of normal geo with vertical land strip inserted

%% To do
% 1) Change met_em initialization to hybrid level-based and change
% 	- 1st level is surface value for 4d vars. CAN CIRCUMVENT BY USING WPS 3.6 AND HISTORICAL 2001 RAP
% 	- SPECHMD (THIS IS QV)
% 	- Soil layers (turn on LSM) - MAKE SURE SST AND SURFACE TEMP MATCH LES
% 2) Standard atmospheric profile above LES
% 3) Make netCDF write block into a function so you can call it in a loop for other files -- note you have to run WPS for all 24 hours on June 7, 2001 to get the radiative forcing correct

%% Write sounding to met_em file
writeIdealMetEm(readFile, writeFile, LESps);
% 
% writeIdealMetEm(readFile2, writeFile2, LESps);
% 
% 
% writeIdealMetEm(readFile3, writeFile3, LESps);
% 
% writeIdealMetEm(readFile4, writeFile4, LESps);
% 
% writeIdealMetEm(readFile5, writeFile5, LESps);

% writeFile1 = writeFile2;
% a = ncread(writeFile1, 'Times');
% a(13) = char('9');
% ncwrite(writeFile1, 'Times', a)

%% CBL
% % Idealized profile from Nieuwstadt et al. (1992) and Soares et al. (2004), as used in the EDMF paper of Witek et al. (2010) [including Teixeira]
% 
% theta_sfc = 300.1; % ??? Based on their figure at initialization
% 
% % z <= 1350
% theta = 300; % [K]
% dqdz = -3.7e-4; % kg / kg / km
% 
% % z > 1350
% dthetadz = 2; % K / km
% dqdz = -9.4e-4; % kg / kg / km
% 
% qsfc = 0.005 * 1000; % g / g converted to kg / kg
% 
% qfx = 2.5e-5; % m kg / kg / s
% hfx = 0.03; % 0.06, 0.09, 0.12 [K m / s]
% 
% u0 = 0.01;
% v0 = 0;
%% Write netCDF file syntax
% 
% % Need to change num_metgrid_levels in info.Dimensions
% % Need to change vertical levels in 3D variables
% % Is it better to use XLONG, XLAT, Angles from AWIP32, then fields from AWP130?
% 
% info = ncinfo([path file1]);
% 
% variables = {'Times', 'PRES', 'GHT', 'SKINTEMP', 'PSFC', 'RH', 'VV', 'UU', 'TT', 'PMSL'};
% 
% 
% 	for var = 1:length(info.Variables)
% 		ncwriteschema(writeFile, info.Variables(var));
% 		matchVar = strcmp( info.Variables(var).Name, variables );
% 		if any( matchVar ) && (flagSkip(matchVar) == 0) % Match between variable and modified variable
% 			ncwrite(writeFile, info.Variables(var).Name, sub.(variables{matchVar})); % Copy changed variable into new netCDF file
% 		else
% 			try
% 				ncwrite(writeFile, info.Variables(var).Name, ncread(fname, info.Variables(var).Name)); % Copy unchanged variable into new netCDF file
% 			catch err
% 			end
% 		end
% 	end
