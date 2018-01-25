% LES data processing code
% Author: Handa Yang - University of California, San Diego
clear all

%% Setup
% Read config file
defaultConfig; divider(1:127) = '=';
% Deal with directories
if ~exist(conf.inputDir, 'dir'), error('Input directory does not exist'), end, if ~exist(conf.outputDir, 'dir'), mkdir(conf.outputDir), end
% Making sure I don't use the wrong directories
fprintf('\n%s\nRun name: %s\nInput directory: %s\nOutput directory: %s\nAveraging interval: %d s (%.1f hr)\n%s\n\n', divider, conf.runName, conf.inputDir, conf.outputDir, conf.avgInterval, conf.avgInterval./3600, divider)
% Plotting setup
lw = 2; fs.axes = 16; fs.legend = 14; style = {'k', 'm', 'b', 'r', 'g', 'c', 'y', 'k--', 'm--', 'b--', 'r--', 'g--', 'c--', 'y--', 'k-.', 'm-.', 'b-.', 'r-.', 'g-.', 'c-.', 'y-.'};
% File saving setup
eps.format = 'epsc2'; png.format = 'png'; fig.format = 'fig';
%% Load raw field variables
% runName.nc: Various 3D fields, sampled at times in 'time'
% Default variable dimensions: u(z, x, y, time)

%% Load profile statistics (1D vars)
% runName.ps.nc: Profile statistics
varNames = [{'time'}, {'zt'}, {'zm'}, {'t'}, {'tot_uw'}, {'sfs_uw'}, {'tot_vw'}, {'sfs_vw'}, {'tot_tw'}, {'sfs_tw'}, {'tot_ww'}, {'sfs_ww'}, {'tot_qw'}, {'sfs_qw'}, {'tot_lw'}, {'sfs_tke'}, {'u_2'}, {'v_2'}, {'w_2'}, {'q'}, {'l'}, {'p'}, {'prd_uw'}, {'u'}, {'v'}, {'boy_prd'}, {'sfs_boy'}, {'sflx'}, {'rflx'}, {'lflxu'}, {'lflxd'}, {'sflxu'}, {'sflxd'}, {'kh'}]; % q in [g/kg]
% zt: t, u_2, v_2, q
% zm: tot_uw, tot_vw, tot_tw, w_2
for idx = 1:length(varNames)
	ps.(varNames{idx}) = ncread([conf.inputDir '/' conf.runName '.ps.nc'], varNames{idx});
end

ps.dt = ps.time(2) - ps.time(1); % Define delta time for temporal averaging purposes
ps.intervalIdx = conf.avgInterval ./ ps.dt;

% Average temporally
for idx = 4:length(varNames) % Omit time, zt, zm from averaging (first three variable names)
	counter = 1;
	for idx2 = 1 : ps.intervalIdx : (size(ps.time, 1) - ps.intervalIdx)
		ps.avg.(varNames{idx}){counter} = mean( ps.(varNames{idx})( :, idx2 : (idx2+ps.intervalIdx) ),  2);
		counter = counter + 1;
	end
end

counter = 1;
for idx = 1 : ps.intervalIdx : (size(ps.time, 1) - ps.intervalIdx)
	ps.segment.stime(counter) = ps.time(idx);
	ps.segment.etime(counter) = ps.time(idx+ps.intervalIdx);
	counter = counter + 1;
end

%% Load time series statistics (Domain averaged vars)
% runName.ts.nc: Time series statistics
% BL depth:
% z_i can be defined in 3 different ways
% * Maximum gradient in Theta_l
% * Maximum variance in Theta_l
% * Minimum buoyancy flux
% NOTE TO SELF: It is possible these are zi1_bar, zi2_bar, and zi3_bar, respectively, as found in the TS file

varNames = [{'time'}, {'zi_bar'}, {'zi1_bar'}, {'zi2_bar'}, {'zi3_bar'}, {'tkeint'}, {'tsrf'}, {'ustar'}, {'vtke'}, {'zcmn'}, {'zbmn'}, {'lwp_bar'}, {'wvp_bar'}, {'shf_bar'}, {'lhf_bar'}, {'cfrac'}, {'sfcbflx'}, {'sflxut'}, {'sflxdt'}, {'sflxus'}, {'sflxds'}, {'lflxut'}, {'lflxdt'}, {'lflxus'}, {'lflxds'}];
for idx = 1:length(varNames)
	ts.(varNames{idx}) = ncread([conf.inputDir '/' conf.runName '.ts.nc'], varNames{idx});
end

ts.dt = ts.time(2) - ts.time(1); % Define delta time for temporal averaging purposes
ts.intervalIdx = conf.avgInterval ./ ts.dt;

% Average temporally
for idx = 2:length(varNames) % Omit time from averaging (first variable name)
	counter = 1;
	for idx2 = 1 : ts.intervalIdx : (size(ts.time, 1) - ts.intervalIdx)
		ts.avg.(varNames{idx}){counter} = mean( ts.(varNames{idx})( idx2 : (idx2+ts.intervalIdx) ),  1);
		counter = counter + 1;
	end
end

counter = 1;
for idx = 1 : ts.intervalIdx : (size(ts.time, 1) - ts.intervalIdx)
	ts.segment.stime(counter) = ts.time(idx);
	ts.segment.etime(counter) = ts.time(idx+ts.intervalIdx);
	counter = counter + 1;
end

% Note to self: u* in TS output; w* = [ (g/T0) * zi * <w * Theta0> ] ^ ⅓, where <w * Theta0> is surface buoyancy flux.
cp_dry = 1006; % [ J / (kg K) ] Specific heat of air at constant pressure
R = 287; % [ J / (kg K) ] Ideal gas constant for air
lambda = 2.5 * 10^6; % [J / kg] Latent heat of vaporization of water
for idx = 1:length(ts.avg.ustar)
	rho = ps.avg.p{idx}(1) ./ (R .* ts.avg.tsrf{idx} );
	ts.avg.wstar{idx} = ( (9.81 ./ ts.avg.tsrf{idx}) .* ts.avg.(conf.zScale){idx} .* ( ps.avg.tot_tw{idx}(1) ./ (cp_dry * rho) ) ) .^ (1/3);
	ts.avg.wstar_sfc{idx} = ( (9.81 .* ts.avg.(conf.zScale){idx} ./ ts.avg.tsrf{idx} ) .* ts.avg.sfcbflx{idx} ) .^ (1/3);
end

% !!! TEMPORARY !!! DELETE THIS
% THIS IS TO SET VALUES TO ABSOLUTE VALUES TO COMPARE WITH WRF-SCM
for idx = 1:length(ts.avg.zi_bar)
	ts.avg.(conf.uScale){idx} = 1;
end

%% Load FIELDS
% 8 processors. 2x4 grid spacing. Large dimension along x-dir, two partitions. Small dimension along y-dir, four partitions.
% Dimensions: z, x, y, time
% Code this way to conserve memory
if ~exist('field', 'var') && isfield(conf, 'loadFields') && (conf.loadFields == 1)
	tmp.fileExt = {'00000000', '00000001', '00000002', '00000003', '00010000', '00010001', '00010002', '00010003'};
	
	% Start, end indices for last hour of data
	field.time = ncread([conf.inputDir '/' conf.runName '.00000000.nc'], 'time');
	field.time_hr = field.time ./ 3600;
	
	field.endIdx = length(field.time_hr);
	field.startIdx = find( field.time_hr == ( (field.time_hr(field.endIdx) - 1) ) );
	
	% Find z/zi = 0.25
	[~, zIdx] = min( abs(ps.zm./ts.avg.(conf.zScale){end} - 0.25) );
	
	% Top side
	for idx = 1:4
		tmp.raw_u = ncread([conf.inputDir '/' conf.runName '.' tmp.fileExt{idx} '.nc'], 'u');
		tmp.raw_w = ncread([conf.inputDir '/' conf.runName '.' tmp.fileExt{idx} '.nc'], 'w');
		tmp.xm = ncread([conf.inputDir '/' conf.runName '.' tmp.fileExt{idx} '.nc'], 'xm');
		tmp.ym = ncread([conf.inputDir '/' conf.runName '.' tmp.fileExt{idx} '.nc'], 'ym');
		tmp.chunk_u(:,:,:) = tmp.raw_u(zIdx, :, :, field.startIdx:field.endIdx);
		tmp.chunk_w(:,:,:) = tmp.raw_w(zIdx, :, :, field.startIdx:field.endIdx);
		if idx == 1
			tmpT.u = tmp.chunk_u;
			tmpT.w = tmp.chunk_w;
			tmpT.xm = tmp.xm;
			tmpT.ym = tmp.ym;
		else
			tmpT.u = horzcat( tmpT.u, tmp.chunk_u );
			tmpT.w = horzcat( tmpT.w, tmp.chunk_w );
			tmpT.xm = horzcat( tmpT.xm, tmp.xm);
			tmpT.ym = horzcat( tmpT.ym, tmp.ym);
		end
		fprintf('Finished loop %d of 8.', idx)
	end
	% Bottom side
	for idx = 5:8
		tmp.raw_u = ncread([conf.inputDir '/' conf.runName '.' tmp.fileExt{idx} '.nc'], 'u');
		tmp.raw_w = ncread([conf.inputDir '/' conf.runName '.' tmp.fileExt{idx} '.nc'], 'w');
		tmp.xm = ncread([conf.inputDir '/' conf.runName '.' tmp.fileExt{idx} '.nc'], 'xm');
		tmp.ym = ncread([conf.inputDir '/' conf.runName '.' tmp.fileExt{idx} '.nc'], 'ym');
		tmp.chunk_u(:,:,:) = tmp.raw_u(zIdx, :, :, field.startIdx:field.endIdx);
		tmp.chunk_w(:,:,:) = tmp.raw_w(zIdx, :, :, field.startIdx:field.endIdx);
		if idx == 5
			tmpB.u = tmp.chunk_u;
			tmpB.w = tmp.chunk_w;
			tmpB.xm = tmp.xm;
			tmpB.ym = tmp.ym;
		else
			tmpB.u = horzcat( tmpB.u, tmp.chunk_u );
			tmpB.w = horzcat( tmpB.w, tmp.chunk_w );
			tmpB.xm = horzcat( tmpB.xm, tmp.xm);
			tmpB.ym = horzcat( tmpB.ym, tmp.ym);
		end
		fprintf('Finished loop %d of 8.', idx)
	end
	
	% Stitch data for LAST HOUR
	field.u = vertcat(tmpT.u, tmpB.u);
	field.w = vertcat(tmpT.w, tmpB.w);
	field.xm = vertcat(tmpT.xm, tmpB.xm);
	field.ym = vertcat(tmpT.ym, tmpB.ym);
	field.note = 'First dim = x, second dim = y. For plotting, transpose and set ydir to normal.';
	field.note2 = 'ym is screwed up. Fix it later--not a priority.';
end

clear tmp tmpT tmpB

%% Main loop

% In neutral conditions, velocity and length scales are u* and zi (there is not necessarily an inversion so zi could be just BL depth).
% In dry convective conditions they are w* and zi.
% And sensible heatfluxes should be normalized by the surface values.
% In wet convective there may be other variables that make sense.
% 
% (<u’w’>^2 + <v’w’>^2)^¼ / u*_sfc

action = 1;

while action == 1
	
	% Options menu
	disp(' '), disp('===== Basic variables ====='), disp(' ')
	disp('[1] Theta_l vs. height')
	disp('[2] q & l vs. height')
	disp('[7] TKE vs. time')
	disp('[8] TKE vs. height')
	
	disp(' '), disp('===== Vertical fluxes ====='), disp(' ')
	disp('[3] Shear stress <u''w''>./w*^2 vs. height')
	disp('[4] Total vertical flux of theta <w''theta''> vs. height')
	disp('[18] Total vertical flux of w <w''w''> vs. height')
	disp('[19] Total vertical flux of q <q''w''> vs. height (WIP)')
	disp('[20] Total vertical flux of l <l''w''> vs. height (WIP)')
	
	disp(' '), disp('===== Cloud properties ====='), disp(' ')
	disp('[11] Cloud fraction vs. time')
	disp('[15] Cloud base and top vs. time')
	disp('[16] LWP vs. time')
	
	disp(' '), disp('===== Radiative fluxes ====='), disp(' ')
	disp('[13] TOA SW vs. time')
	disp('[14] TOA LW vs. time')
	disp('[22] Radiative flux vs. height')
	disp('[23] LW flux vs. height')
	disp('[24] SW flux vs. height')
	
	disp(' '), disp('===== Velocity ====='), disp(' ')
	disp('[9] Velocity variance vs. time')
	disp('[10] Velocity variance vs. height')

	disp(' '), disp('===== Surface fluxes and buoyancy ====='), disp(' ')
	disp('[12] Surface buoyancy flux vs. time')
	disp('[17] Surface heat flux vs. time')
	disp('[21] Buoyancy production of resolved TKE vs. height (WIP)')
	
	disp(' '), disp('===== Misc. ====='), disp(' ')
	disp('[99] Plot shear stress (vector averaged) <u''w''>, <v''w''>')
	
	disp('[0] Ensemble')
	
	% Prompt input
	disp(' ')
	selection = input('Input menu item: ');
	
	if isempty(selection), continue, end
	
	% Deal with selected option
	switch(selection)
		case 1 % Theta (PS) vs. height
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(1); plotTheta(ps.avg.t, ps.segment, ps.zt, ts.avg.(conf.zScale), plotSegments, lw, fs, style)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_theta' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_theta' '.png'], png); hgexport(h, [conf.outputDir '/fig_' conf.runName '_theta' '.fig'], fig);
			
		case 2 % q (PS) & l (PS) [Water vapor mixing ratio, liquid water mixing ratio] vs. height
			
			%%%%% NOTE TO SELF: ADD L TO THIS PLOT %%%%%
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(2); plotWater(ps, ts, conf, plotSegments, lw, fs, style)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_water' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_water' '.png'], png);
		
		case 3 % Shear stress (PS) (Sum of resolved + SGS) vs. height
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(3); plotShearStress(ps, ts, ps.segment, ps.zm, ts.avg.(conf.zScale), plotSegments, conf.uScale, lw, fs, style)
			
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_shearStress' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_shearStress' '.png'], png);
			
		case 4 % Heat flux (PS) (Sum of resolved + SGS) vs. height
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(4); plotHeatFlux(ps, ps.segment, ps.zm, ts.avg.(conf.zScale), plotSegments, lw, fs, style)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_heatFlux' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_heatFlux' '.png'], png);
            
        case 5 % Km
            
            showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
            h = newFig(5); plotTheta(ps.avg.kh, ps.segment, ps.zt, ts.avg.(conf.zScale), plotSegments, lw, fs, style)
            xlabel('K_m [m^2/s]')
            
            hgexport(h, [conf.outputDir '/eps_' conf.runName '_Kh' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_Kh' '.png'], png); hgexport(h, [conf.outputDir '/fig_' conf.runName '_Kh' '.fig'], fig);
            
		case 7 % TKE vs. time
			
			h = newFig(1); plotTKE(ts.time, ts.tkeint, lw, fs)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_TKE' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_TKE' '.png'], png);
			
		case 8 % TKE vs. height
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(2); plotTKEprofile(ps, ts, conf, plotSegments, lw, fs, style)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_TKE' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_TKE' '.png'], png);
			
		case 9 % Velocity variance vs. time (PS)
			
			h = newFig(3); plotVar(ps, lw, fs)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_velVar' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_velVar' '.png'], png);
			
		case 10 % Velocity variance vs. height (PS)
			
			h(1) = newFig(4); h(2) = newFig(5); h(3) = newFig(6); plotSegments = input('Input which times to plot: ');
			plotVarProfiles(ps, ts, conf, plotSegments, lw, fs, style, h)
			
			hgexport(h(1), [conf.outputDir '/eps_' conf.runName '_velVar_u' '.eps'], eps); hgexport(h(1), [conf.outputDir '/png_' conf.runName '_velVar_u' '.png'], png);
			hgexport(h(2), [conf.outputDir '/eps_' conf.runName '_velVar_v' '.eps'], eps); hgexport(h(2), [conf.outputDir '/png_' conf.runName '_velVar_v' '.png'], png);
			hgexport(h(3), [conf.outputDir '/eps_' conf.runName '_velVar_w' '.eps'], eps); hgexport(h(3), [conf.outputDir '/png_' conf.runName '_velVar_w' '.png'], png);
			
		case 11 % Cloud fraction vs. time
			
			h = newFig(1); plotCloudFrac(ts, lw, fs)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_CloudFrac' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_TKE' '.png'], png);
			
		case 12 % Surface buoyancy flux vs. time
			
			h = newFig(1); plotSurfBuoyancy(ts, lw, fs)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_SurfBuoyancy' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_TKE' '.png'], png);
			
		case 13 % Top of atmosphere shortwave radiation vs. time
			
			h = newFig(1); plotSW(ts, lw, fs)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_TOA_SW' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_TKE' '.png'], png);
			
		case 14 % Top of atmosphere longwave radiation vs. time
			
			h = newFig(1); plotLW(ts, lw, fs)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_TOA_LW' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_TKE' '.png'], png);
			
		case 15 % Cloud base and top height vs. time
			
			h = newFig(1); plotCloudHeight(ts, lw, fs)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_CloudHeight' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_TKE' '.png'], png);
			
		case 16 % Liquid water path vs. time
			
			h = newFig(1); plotLWP(ts, lw, fs)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_LWP' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_TKE' '.png'], png);
			
		case 17 % Heat flux vs. time (sensible, latent)
			
			h = newFig(1); plotSurfHeatFlux(ts, lw, fs)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_SurfHeatFlux' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_TKE' '.png'], png);
			
		case 18 % Total vertical flux of w vs. height
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(3); plotFluxW(ps, ts, conf, plotSegments, lw, fs, style)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_fluxW' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_shearStress' '.png'], png);
			
		case 19 % Total vertical flux of q vs. height
			
			%% WIP %%
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(3); % plotFluxQ(ps, ts, conf, plotSegments, lw, fs, style);
            plotFluxQ(ps, ps.segment, ps.zm, ts.avg.(conf.zScale), plotSegments, lw, fs, style);
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_fluxQ' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_shearStress' '.png'], png);
			
		case 20 % Total vertical flux of l (not including sfs?) vs. height
			
			%% WIP %%
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(3); plotFluxL(ps, ts, conf, plotSegments, lw, fs, style)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_fluxL' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_shearStress' '.png'], png);
			
		case 21 % Buoyancy production of resolved TKE vs. height
			
			%% WIP
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(3); plotBuoyProd(ps, ts, conf, plotSegments, lw, fs, style)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_buoyProd' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_shearStress' '.png'], png);
			
		case 22 % Radiative flux vs. height (total, shortwave, longwave [subtract short from total])
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(3); plotRadFlux(ps, ts, conf, plotSegments, lw, fs, style)
						
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_radFlux' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_shearStress' '.png'], png);
			
		case 23 % Longwave radiative flux up and down vs. height
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(3); plotRadLW(ps, ts, conf, plotSegments, lw, fs, style)
						
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_radLW' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_shearStress' '.png'], png);
			
		case 24 % Shortwave radiative flux up and down vs. height
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(3); plotRadSW(ps, ts, conf, plotSegments, lw, fs, style)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_radSW' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_shearStress' '.png'], png);
			
		case 25 % Mean x-direction velocity vs. height
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(3); plotMeanU(ps, ts, conf, plotSegments, lw, fs, style)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_meanU' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_shearStress' '.png'], png);
			
		case 99 % Shear stress (PS) (Sum of resolved + SGS)
			
			showSegments(ps.segment, conf.avgInterval, divider); plotSegments = input('Input which times to plot: ');
			h = newFig(12); plotShearStress(ps, ts, ps.segment, ps.zm, ts.avg.(conf.zScale), plotSegments, conf.uScale, lw, fs, style, 1)
			
			hgexport(h, [conf.outputDir '/eps_' conf.runName '_shearStress' '.eps'], eps); hgexport(h, [conf.outputDir '/png_' conf.runName '_shearStress_vecAvg' '.png'], png);
			
		case 0 % Ensemble of 1-6
			
			h = newFig(99);
			
			hgexport(h(1), [conf.outputDir '/eps_' conf.runName '_ensemble' '.eps'], eps); hgexport(h(1), [conf.outputDir '/png_' conf.runName '_ensemble' '.png'], png);
			
		otherwise % Catch to prevent infinite loop
			
			action = 0;
			
	end
	
end

%% Testing stuff
figure
for z = 1:24
	plot(ps.avg.tot_vw{z}, ps.zm)
	xlim([-0.1 0.1]), ylim([0 1200])
	pause(0.3)
end

%% GIBBS
for z = 1:61, imagesc(field.u(:,:,z)', [-6 6]), set(gca, 'ydir', 'normal'), colormap(redblue), title(sprintf('%d', z)), pause, end
for z = 1:61, imagesc(field.w(:,:,z)', [-6 6]), set(gca, 'ydir', 'normal'), colormap(redblue), title(sprintf('%d', z)), pause, end
% Mean u
figure; plot(ps.avg.u{12}./ts.avg.wstar{12}, ps.zm./ts.avg.zi1_bar{12}, 'k', 'linewidth', 2), set(gca,'fontsize', fs.axes)

% 1D spectral
figure; uSpectra_k(field, ts, conf, lw, fs); figure; uSpectra_k2(field, ts, conf, lw, fs); figure; wSpectra_k(field, ts, conf, lw, fs); figure; wSpectra_k2(field, ts, conf, lw, fs);

% 2D spectral

twoDimSpectra_u(field, ts, conf, fs), twoDimSpectra_w(field, ts, conf, fs)

% Do video!
slice.xz = ncread([conf.inputDir '/' conf.runName '.out.xz.0000.0002.nc'], 't');
slice.xz(:,129:256,:) = ncread([conf.inputDir '/' conf.runName '.out.xz.0001.0002.nc'], 't');
%%
h = figure('position', [200 200 1000 220]);
for idx = 1:size(slice.xz,3)
	imagesc(field.xm(:,1), ps.zm, slice.xz(:,:,idx), [300 310])
	title(sprintf('t = %.1f hr', (idx-1)*60/3600))
	set(gca, 'ydir', 'normal', 'fontsize', fs.axes)
	colorbar
	set(gca, 'ydir', 'normal', 'fontsize', fs.axes)
	xlabel('Distance [m]', 'fontsize', fs.axes)
	ylabel('Height [m]', 'fontsize', fs.axes)
% 	pause(0.005)
% 	pause
% 	hgexport(h, ['shear_40/video/' sprintf('%04d', idx) '.png'], png)
	hgexport(h, ['free_40/video2/' sprintf('%04d', idx) '.png'], png)
end
