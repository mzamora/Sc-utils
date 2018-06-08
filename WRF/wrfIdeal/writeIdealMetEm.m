%% Make shell of netCDF file

function [out] = writeIdealMetEm(readFile1, writeFile1, LESps)

% Constants
defineConstants

%% Get rid of old file if present
delete(writeFile1)

%% Load LES variables
lesVar = {'time', 'p', 'RH', 'q', 'l', 'zm', 'v', 'u', 't'}; % t = theta_l, so needs processing to put it into theta, and qv = (q - l) [q = qt]

for n = 1:length(lesVar)
	les.(lesVar{n}) = ncread(LESps, lesVar{n});
end

%% LES Postprocessing

% DYCOMS
les.RH = les.RH.*100; % Put in % instead of decimal
les.theta = les.t + Lv / cpAir .* les.l ./ 1000; % les.t is theta_l
les.temperature = tFromTheta( les.theta, les.p );
les.qv = (les.q - les.l) ./ 1000; % les.q is qt
les.l = les.l ./ 1000; % Mixing ratios in LES are in g/kg
les.zmOrig = les.zm;

% Add standard atmosphere above LES

%%%%%% THIS IS OFFSET BY 1 LEVEL -- Z AT 149 TO TOA
%%%%%% 

% INVTOPLEVEL = 149; % DYCOMS RF01

[dT_max, numofInv, ziTop, ziBase, ~, ~]=TMP_Inversion_Strength_Cal_V1(les.temperature(:,1), les.zm./1000, 0);

INVTOPLEVEL = findMatch(les.zm, ziBase*1000)+1; % 1 point above zi Base for LES high resolution profiles

%%%
figure; plot(les.temperature(:,1), les.zm)
hold on
%%%

zstd = (les.zm(INVTOPLEVEL,1)+200):200:21000;

% Determine offset from standard atmosphere
[~,~,T_originalStandard,~,~,~] = atmos(les.zm(INVTOPLEVEL), 'altType', 'geometric', 'tOffset',0);
tOffset = les.temperature(INVTOPLEVEL) - T_originalStandard;

[~,~,Tstd,Pstd,~,~] = atmos(zstd, 'altType', 'geometric', 'tOffset', tOffset); % 21000 meters is below 50 hPa [DYCOMS]
% [~,~,Tstd,Pstd,~,~] = atmos(zstd, 'altType', 'geometric', 'tOffset', 5); % 21000 meters is below 50 hPa [SIEBESMA]
% [~,~,Tstd,Pstd,~,~] = atmos(zstd, 'altType', 'geometric', 'tOffset', 7.7); % 21000 meters is below 50 hPa [WITEK]

les.zm = [les.zm(1:INVTOPLEVEL); zstd'];
les.RH = [les.RH(1:INVTOPLEVEL,:); repmat(8, length(zstd), size(les.p,2))];
les.temperature = [les.temperature(1:INVTOPLEVEL,:); repmat(Tstd', 1, size(les.temperature,2))]; % Need to repmat this to fill in bottom of temperature matrix
les.t = les.temperature; % Replace t (theta_l) with temperature
les.p = [les.p(1:INVTOPLEVEL,:); repmat(Pstd', 1, size(les.p,2))];

%%%
plot(les.temperature(:,1), les.zm, 'r--')
%%%

% DYCOMS: Make logarithmic velocity profiles
% les.u = []; les.v = []; % Simpler to delete u and v
% Exponent can be 1/7 in neutrally stable conditions over open land. 0.11 is more appropriate over open ocean.
% % DYCOMS RF01
% les.u = 6 .* (les.zm ./ 500) .^ (0.11); % Use 0.11 as exponent for open ocean. 6 m/s observed in DYCOMS II RF01. Assume at 500 m (middleish of PBL).
% les.v = -4.25 .* (les.zm ./ 500) .^ (0.11); % Use 0.11 as exponent for open ocean. -4.25 m/s observed in DYCOMS II RF01. Assume at 500 m (middleish of PBL).
% % les.u = les.u(:,1); les.v = les.v(:,1);
% les.u(INVTOPLEVEL:end) = 7; % U_g
% les.v(INVTOPLEVEL:end) = -5.5; % V_g

% CGILS S12
% les.u = 1 .* (les.zm ./ 500) .^ (0.11); % Use 0.11 as exponent for open ocean. 6 m/s observed in DYCOMS II RF01. Assume at 500 m (middleish of PBL).
% les.v = -10 .* (les.zm ./ 500) .^ (0.11); % Use 0.11 as exponent for open ocean. -4.25 m/s observed in DYCOMS II RF01. Assume at 500 m (middleish of PBL).
% les.u(INVTOPLEVEL:end) = 1; % U_g
% les.v(INVTOPLEVEL:end) = -10; % V_g
% Try straight interpolation
origUTop = les.u(end,1);
origVTop = les.v(end,1);
les.u = interp1(les.zmOrig, les.u(:,1), les.zm, 'linear', origUTop);
les.v = interp1(les.zmOrig, les.v(:,1), les.zm, 'linear', origVTop);

% CBL
% les.u = 6 .* (les.zm ./ 675) .^ (1./7); % try dycoms?
% les.v = -4.25 .* (les.zm ./ 675) .^ (1./7); % Try with dycoms velocity
% les.u = zeros(size(les.zm)); %0.01 .* (les.zm ./ 675) .^ (1./7); % Use 0.11 as exponent for open ocean. 6 m/s observed in DYCOMS II RF01. Assume at 500 m (middleish of PBL).
% les.v = zeros(size(les.zm)); % Use 0.11 as exponent for open ocean. -4.25 m/s observed in DYCOMS II RF01. Assume at 500 m (middleish of PBL).
% les.u(INVTOPLEVEL:end) = 0; % U_g
% les.v(INVTOPLEVEL:end) = 0; % V_g

les.u = repmat(les.u, [1 size(les.zm,2)]);
les.v = repmat(les.v, [1 size(les.zm,2)]);

% Fill above atmosphere stuff
for n = size(les.qv, 1):size(les.RH,1)
	les.qv(n, :) = les.RH(n, :) .* qsat(les.t(n, :), les.p(n, :)) ./ 100;
% 	les.u(n, :) = les.u(n-1, :); % Put previous values back in u and v
% 	les.v(n, :) = les.v(n-1, :);
	les.l(n, :) = 0;
end


% les.qv = [les.qv; repmat((linspace(0.015, 2e-6, length(zstd)))', 1, size(les.qv,2))];

% Make matrices for QG (graupel), QS (snow), QI (ice), QR (rain)
les.QG = zeros(size(les.t));
les.QS = zeros(size(les.t));
les.QI = zeros(size(les.t));
les.QR = zeros(size(les.t));

% % % % % % % % % % % % Quick hack to fix top level
% % % % % % % % % % % % Long-term solution is to implement standard atmosphere for top of model
% % % % % % % % % % % les.p(end,1) = 28344.1;
% % % % % % % % % % % les.t(end,1) = tFromTheta( 337.845, les.p(216,1));
% % % % % % % % % % % les.RH(end,1) = 50;
% % % % % % % % % % % les.zm(end,1) = 98214.3;

% Reorder LES data such that it matches hybrid grid RAP
lesVar = {'time', 'p', 'RH', 'QG', 'QS', 'QI', 'QR', 'l', 'p', 'zm', 'v', 'u', 'qv', 't'}; % t = theta_l, so needs processing to put it into theta qv = (q - l)
%%%%%%%% NOTE IT IS IMPORTANT THAT THE ORDER IS MAINTAINED %%%%%%%%%%%%%

% Vertical levels in LES data
numLESvert = size(les.p, 1);

% Read netCDF schema info
info = ncinfo(readFile1);
edited = info;

% Replace num_metgrid_levels in dimension definitions with LES resolution
edited.Dimensions(5).Length = numLESvert;
edited.Attributes(5).Value = int32(numLESvert);

count = 0;
for n = 2:length(info.Variables) % All variables excluding Times
	% Replace num_metgrid_levels with LES vertical levels -- first search for match in Dimensions name
	for m = 1:length(info.Variables(n).Dimensions)
		if strcmpi(info.Variables(n).Dimensions(m).Name, 'num_metgrid_levels')
			% Found vertical variable--replace with LES dimensions
			edited.Variables(n).Dimensions(m).Length = numLESvert;
			edited.Variables(n).Size(m) = numLESvert;
			% Print found variables and store in a cell array
			fprintf('Replaced num_metgrid_levels for %s.\n', info.Variables(n).Name)
			count = count + 1;
			profileVar{count} = info.Variables(n).Name;
			profileVarId{count} = n;
		end
	end
end

% Corresponding LES variables -- 1 is always time, so skip it.

for n = 1:length(profileVarId)
	for i = 1:edited.Variables(profileVarId{n}).Size(1)
		for j = 1:edited.Variables(profileVarId{n}).Size(2)
			writeVar{n}(i, j, :) = les.(lesVar{n+1})(:,1); % 1 is always time, so skip it.
		end
	end
end

% variables = {'Times', 'PRES', 'GHT', 'SKINTEMP', 'PSFC', 'RH', 'VV', 'UU', 'TT', 'PMSL'};

% Create structure of new netCDF file
ncwriteschema(writeFile1, edited);

% Write variables to new netCDF file
count = 0;
for n = 1:length(info.Variables)
	if any(n == [profileVarId{:}])
		count = count+1;
		ncwrite(writeFile1, info.Variables(n).Name, writeVar{count});
		fprintf('Writing modified %s in %s.\n', info.Variables(n).Name, writeFile1)
	else
		ncwrite(writeFile1, info.Variables(n).Name, ncread(readFile1, info.Variables(n).Name));
		fprintf('Writing original %s in %s.\n', info.Variables(n).Name, writeFile1)
	end
end

% Write 2D matrices PSFC and PMSL
spatialVar = {'PSFC', 'PMSL'};
for n = 1:length(spatialVar)
	ncwrite(writeFile1, spatialVar{n}, repmat(les.p(1,1),edited.Dimensions(3).Length,edited.Dimensions(4).Length)); % Homogeneous PSFC and PMSL taken from 1st LES model level
% 	ncwrite(writeFile1, spatialVar{n}, repmat(100000,edited.Dimensions(3).Length,edited.Dimensions(4).Length)); % Homogeneous PSFC and PMSL taken from 1st LES model level
	fprintf('Writing modified %s in %s.\n', spatialVar{n}, writeFile1)
end

% CBL
% ncwrite(writeFile1, 'SKINTEMP', repmat(300,edited.Dimensions(3).Length,edited.Dimensions(4).Length)); % Homogeneous PSFC and PMSL taken from 1st LES model level
% ncwrite(writeFile1, 'LU_INDEX', repmat(6,edited.Dimensions(3).Length,edited.Dimensions(4).Length)); % Homogeneous PSFC and PMSL taken from 1st LES model level

%% CBL REPLACE SURFACE PARAMETERS
i = 25; j = 25; % Land pixel that we're going to copy everything from
% Search 3D
surface3D = {'HGT_M', 'LANDMASK', 'LANDSEA', 'OA1', 'OA2', 'OA3', 'OA4', 'OL1', 'OL2', 'OL3', 'OL4', 'SCB_DOM', 'SCT_DOM', 'SLOPECAT', 'SOILM000', 'SOILM160', 'SOILM300', ...
	'SOILT000', 'SOILT160', 'SOILT300', 'SOILTEMP', 'VAR', 'VAR_SSO', 'LU_INDEX'};
for n = 1:length(info.Variables)
	if any(strcmpi(info.Variables(n).Name, surface3D))
		tempColumn = ncread(readFile1, info.Variables(n).Name);
		tempColumn = tempColumn(i,j);
		replacementVar = repmat(tempColumn,[edited.Dimensions(3).Length,edited.Dimensions(4).Length]);
		if strcmpi(info.Variables(n).Name, 'HGT_M')
			replacementVar(:) = 0;
		end
		% Write new variable
		ncwrite(writeFile1, info.Variables(n).Name, replacementVar);
		fprintf('Writing modified %s in %s.\n', info.Variables(n).Name, writeFile1)
	end
end
% Search 4D
surface4D = {'ALBEDO12M', 'GREENFRAC', 'LAI12M', 'LANDUSEF', 'SOILCBOT', 'SOILCTOP', 'SOILM', 'SOILT'};
for n = 1:length(info.Variables)
	if any(strcmpi(info.Variables(n).Name, surface4D))
		tempColumn = ncread(readFile1, info.Variables(n).Name);
		tempColumn = tempColumn(i,j,:);
		replacementVar = repmat(tempColumn,[edited.Dimensions(3).Length,edited.Dimensions(4).Length]);
		% Write new variable
		ncwrite(writeFile1, info.Variables(n).Name, replacementVar);
		fprintf('Writing modified %s in %s.\n', info.Variables(n).Name, writeFile1)
	end
end

% DYCOMS RF01
% ncwrite(writeFile1, 'SKINTEMP', repmat(292,edited.Dimensions(3).Length,edited.Dimensions(4).Length)); % Homogeneous PSFC and PMSL taken from 1st LES model level
% ncwrite(writeFile1, 'LU_INDEX', repmat(14,edited.Dimensions(3).Length,edited.Dimensions(4).Length)); % SOILPARM.TBL: 14 = water

% CGILS S12 CTL
ncwrite(writeFile1, 'SKINTEMP', repmat(290.6,edited.Dimensions(3).Length,edited.Dimensions(4).Length)); % Homogeneous PSFC and PMSL taken from 1st LES model level
ncwrite(writeFile1, 'LU_INDEX', repmat(14,edited.Dimensions(3).Length,edited.Dimensions(4).Length)); % SOILPARM.TBL: 14 = water
%% ZERO CORIOLIS PARAMETERS
% ncwrite(writeFile1, 'F', zeros(edited.Dimensions(3).Length, edited.Dimensions(4).Length));
% ncwrite(writeFile1, 'E', zeros(edited.Dimensions(3).Length, edited.Dimensions(4).Length));

fprintf('Done.\n')

out = 1;

%% Snippet to write to LES output sounding
% fid = fopen('soundingLES','w');
% 
% fopen(fid);
% fprintf(fid,' 1018.690   %7.3f   %6.3f    0.000    0.000\n', les.theta(1,1), 1000.*les.qv(1,1));
% for n = 1:length(les.zm)
% 	fprintf(fid,'%9.3f   %7.3f   %6.3f   %6.3f   %6.3f\n', les.zm(n), les.theta(n,1), 1000.*les.qv(n,1), les.u(n,1), les.v(n,1));
% end
% fclose(fid);

end