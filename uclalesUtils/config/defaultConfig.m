% ============================================================================
% LES data processing configuration file
% ============================================================================

% ==== CONFIG PARAMETERS =====================================================

% runName: Name of output file prefix
% inputDir: LES raw data.
% outputDir: Plots, etc.
% avgInterval: Time averaging interval in seconds
% zScale: zi_bar, zi1_bar, zi2_bar, zi3_bar (???, maximum grad(Theta_l), max var(Theta_l), min buoyancy flux)
% uScale: ustar, wstar

%%
%%% %%% %%% %%% %%% %%% %%%

%%% %%%     WRF     %%% %%%

%%% %%% %%% %%% %%% %%% %%%

%%

% conf.runName		=		'wrfles_06102013';
% conf.inputDir		= 		'~/uclales/run/wrfles_06102013/';
% conf.outputDir		=		'~/research/lesData/wrfles_06102013/';
% conf.avgInterval	=		1800;
% conf.zScale			=		'zi3_bar';
% conf.uScale			=		'wstar';

%%
%%% %%% %%% %%% %%% %%% %%%

%%% %%%  Idealized  %%% %%%

%%% %%% %%% %%% %%% %%% %%%

%% UCLALES CBL

conf.runName			=	'dcbl_x';
if ispc
conf.inputDir		=	'C:\Users\Handa\Dropbox\matlab\lesData\dcbl_witek_2010\';
conf.outputDir       =   'C:\Users\Handa\Documents\Work\lesData\dcbl_witek\';
elseif ismac
    conf.inputDir = '~/Documents/MATLAB/cbl/';
    conf.outputDir = '~/Documents/MATLAB/cbl/plots/';
end
conf.avgInterval		=	600;
conf.zScale			=	'zi1_bar';
conf.uScale			=	'wstar_sfc';


%% TO COMPARE AGAINST WITEK 2010 CBL
% conf.runName			=	'dcbl_x';
% if ispc
% conf.inputDir		=	'C:\Users\Handa\Dropbox\matlab\lesData\dcbl_witek_2010\';
% conf.outputDir       =   'C:\Users\Handa\Documents\Work\lesData\dcbl_witek\';
% elseif ismac
%     conf.inputDir = '~/Dropbox/matlab/lesData/dcbl_witek_2010/';
%     conf.outputDir = '~/Dropbox/matlab/lesData/dcbl_witek_2010/plots/';
% end
% conf.avgInterval		=	3600;
% conf.zScale			=	'zi1_bar';
% conf.uScale			=	'wstar_sfc';

%%
%runName				dcbl_x
%inputDir			~/uclales/run/dcbl_gibbs_shear_real_40/
%outputDir			~/research/lesData/dcbl_gibbs_shear_real_40/
%inputDir			/media/Windows7_OS/TRANSFER/dcbl_gibbs_shear_40_HF/
%outputDir           ~/research/lesData/GIBBS_SHR_40/
%inputDir			/media/Windows7_OS/TRANSFER/dcbl_gibbs_free_40/
%outputDir           ~/research/lesData/GIBBS_FREE_40/
%avgInterval			3600
%zScale				zi1_bar
%uScale				wstar

%%
% conf.runName				= 'rf01';
% conf.inputDir			= '~/uclales/run/dycom1_div/';
% conf.outputDir			= '~/research/lesData/dycom1_div/';
% conf.avgInterval			= 3600;
% conf.zScale				= 'zi3_bar';
% conf.uScale				= 'ustar';

%% Dissipating Sc
% conf.runName		=		'cgils_s12_ctl';
% conf.inputDir		=		'~/research/matlab/data/LES/';
% conf.outputDir		=		'~/research/matlab/OUTPUT/';
% conf.avgInterval	=		2400;
% conf.avgInterval	=		1200;
% conf.zScale			=		'zi3_bar';
% conf.uScale			=		'wstar';
% conf.uScale			=		'wstar_sfc';

%% CGILS

% conf.runName		= 'cgils_s12_ctl';
% conf.inputDir		= '/mnt/lab_24tb1/users/mghonima/p21/output/24_hr_unmodified_CGILS_specs_solar_varying_strtime_196';
% conf.outputDir		= '~/research/matlab/OUTPUT/';
% conf.avgInterval    = 2400;
% conf.zScale			= 'zi3_bar';
% conf.uScale			= 'wstar_sfc';

%% Dissipating Sc
% conf.runName		=		'cgils_s12_ctl';
% conf.inputDir		=		'~/research/matlab/data/LES/minus1/';
% conf.outputDir		=		'~/research/matlab/OUTPUT/';
% conf.avgInterval	=		2400;
% conf.avgInterval	=		1200;
% conf.zScale			=		'zi1_bar';
% conf.uScale			=		'wstar';
% conf.uScale			=		'wstar_sfc';

%%
%%% %%% %%% %%% %%% %%% %%%

%%% %%% First tries %%% %%%

%%% %%% %%% %%% %%% %%% %%%

%runName				dcbl_x
%inputDir			~/uclales/run/firstTry_dcbl_coarse/
%outputDir			~/research/lesData/firstTrydcbl_coarse/
%avgInterval			3600
%zScale				zi3_bar
%uScale				wstar

%runName				rf01
%inputDir			~/uclales/run/firstTry_dycom1/
%outputDir			~/research/lesData/firstTrydycom1/
%avgInterval			3600
%zScale				zi3_bar
%uScale				ustar

%runName				astx40
%inputDir			~/uclales/run/astex/
%outputDir			~/research/lesData/astex/
%avgInterval		3600
%zScale				zi3_bar
%uScale				ustar