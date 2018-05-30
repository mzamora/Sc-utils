%% Retrieves select variables for a defined CV within WRF domain. Then computes heat and moisture budgets, as well as cloud base tendency using mixed-layer theory.
% [cv] = wrfMLCV(wrf, startTime, endTime, x, y, streamer)
%
% wrf: struct containing WRF variables (native names)
%
% startTime, endTime: datenums defining start and end times in UTC
%
% x, y: Indices in WRF defining CV region.
%       IMPORTANT: Take 1 extra grid point around CV in order to compute spatial derivatives.
%       Determine by looking at an imagesc of WRF data. x and y coordinates correspond to x, y, on data cursor.
%
% streamer: struct containing streamer output. Leave empty if you wish to use WRF radiative theta tendency.
%           NOTE: WRF radiative tendency is THETA tendency. Strictly speaking, the heat budget is used to calculated THETA_L tendency. The error from this approximation is small, though.
%
% NOTE: Due to a stupid coding error, matrices are one size larger than they should be.
%       I.e., if you gave the CV two 6 element vectors (4x4 domain with 1 grid point border for spatial derivatives), the output is a 5x5 matrix rather than a
%		4x4 matrix. The elements in the 1st row and 1st column of the 5x5 matrix will be 0

function [cv] = wrfMLCVavg(wrf, startTime, endTime, x, y, soundingDir, streamer)

if nargin < 6
	soundingDir = [];
	streamer = [];
elseif nargin < 7
	streamer = [];
end

if ~isempty(soundingDir)
	if ~exist(soundingDir, 'dir'), mkdir(soundingDir), end
	makeSounding = 1;
else
	makeSounding = 0;
end

%% Define constants
Lv = 2.5e6; % [J / kg]
densityWater = 1e3; % [kg / m^3]
% densityAir = 1.2922; % [kg / m^3] (Wikipedia)
densityAir = 1.2; % [kg / m^3] (http://www.engineeringtoolbox.com/air-density-specific-weight-d_600.html, approx. at 20 C [68 F])
cpAir = 1005; % [J / kg K] (http://www.engineeringtoolbox.com/air-properties-d_156.html at 20 C)
R = 287; % Specific gas constant for dry air (also called Rd) [J / kg K]
Rv = 461; % Specific gas constant for water vapor [J / kg K]
g = 9.81; % [m / s^2]
D = 3.75 * 10^-6; % [1 / s]

%% WRF SIMULATION CONSTANTS
h = 8100; % dx = dy = 2.5 km = 2500 m

%% CV variable retrieval
simStart = findMatch( round(wrf.time.*86400), round(startTime.*86400));
simEnd =   findMatch( round(wrf.time.*86400), round(endTime.*86400));
cv.time = wrf.time(simStart:simEnd) - 8/24; % PST
cv.znu = wrf.ZNU;
cv.znw = wrf.ZNW;

% Obtain height at theta levels for points within domain
for t = 1:length(cv.time)
	tWRF = t + simStart - 1; % Corresponding time index in WRF time vector
	for i = 1:length(x)
		for j = 1:length(y)
			for z = 1:size(wrf.U, 3)
				cv.heightstag(i, j, z, t) = ( wrf.PH( x(i), y(j), z, tWRF ) + wrf.PHB( x(i), y(j), z, tWRF ) ) ./ 9.81;
				cv.height(i, j, z, t) = 0.5 .* ( wrf.PHB( x(i), y(j), z, tWRF) + wrf.PH( x(i), y(j), z, tWRF) + wrf.PHB( x(i), y(j), z + 1, tWRF) + wrf.PH( x(i), y(j), z + 1, tWRF) ) ./ 9.81; % MSL
				cv.heightAGL(i, j, z, t) = cv.height(i, j, z, t) - wrf.HGT( x(i), y(j), tWRF ); % AGL
			end
		end
	end
end

% Examine spatially averaged mixed layer profiles within CV
for t = 1:length(cv.time)
	tWRF = t + simStart - 1; % Corresponding time index in WRF time vector
    temp = wrf.ZNW(:, tWRF);
	cv.DNW(:, t) = diff(temp);
	for i = 1:length(x)
		for j = 1:length(y)
			cv.spatial.hgt(i, j, t) = wrf.HGT( x(i), y(j), tWRF); % Ground elevation [m]
			try
			cv.spatial.idxPBLH(i, j, t) = findMatch(squeeze(cv.heightAGL(i, j, :, t)), wrf.PBLH(x(i), y(j), tWRF)); % Index corresponding to PBLH variable [-]
			catch
				cv.spatial.idxPBLH(i, j, t) = 1;
			end
			cv.spatial.PBLH(i, j, t) = wrf.HGT( x(i), y(j), tWRF ) + wrf.PBLH( x(i), y(j), tWRF ); % PBL height determined by PBL scheme [m]
			cv.spatial.shf(i, j, t) = wrf.HFX(x(i), y(j), tWRF); % Sensible heat flux [W / m^2]
			cv.spatial.lhf(i, j, t) = wrf.LH(x(i), y(j), tWRF); % Latent heat flux [W / m^2]
			cv.spatial.psfc(i, j, t) = wrf.PSFC(x(i), y(j), tWRF); % Surface pressure [Pa]
			cv.spatial.RainSFC(i, j, t) = wrf.RAINC( x(i), y(j), tWRF) + wrf.RAINNC( x(i), y(j), tWRF); % Total accumulated rain at surface (Cu scheme + grid-scale). Typically starts from 00:00 UTC and is stored for 24 hours [mm].
			cv.spatial.tsk( i, j, t ) = wrf.TSK( x(i), y(j), tWRF ); % Skin temperature of surface [K]
			cv.spatial.q2( i, j, t ) = wrf.Q2( x(i), y(j), tWRF ); % Skin temperature of surface [K]
			cv.spatial.ql2( i, j, t ) = cv.spatial.q2( i, j, t ) - qsat( wrf.T2( x(i), y(j), tWRF), cv.spatial.psfc(i, j, t) ); % Skin temperature of surface [K]
			cv.spatial.tl2( i, j, t ) = wrf.TH2( x(i), y(j), tWRF ) - Lv/cpAir .* (wrf.TH2( x(i), y(j), tWRF ) ./ wrf.T2( x(i), y(j), tWRF )) .* cv.spatial.ql2(i, j, t); % Theta_l at 2m [K]
% 			cv.spatial.swdn( i, j, t ) = wrf.SSWDN( x(i), y(j), tWRF);
% 			cv.spatial.swup( i, j, t ) = wrf.SSWUP( x(i), y(j), tWRF);
% 			cv.spatial.lwdn( i, j, t ) = wrf.SLWDN( x(i), y(j), tWRF);
% 			cv.spatial.lwup( i, j, t ) = wrf.SLWUP( x(i), y(j), tWRF);
			for z = 1:size(cv.height, 3)
				% Interpolate all variables onto theta levels (grid point centers) [u, v, w]
				cv.spatial.pressure( i, j, z, t ) = wrf.PB( x(i), y(j), z, tWRF); % Atmospheric pressure [Pa]
				cv.spatial.qt( i, j, z, t ) = wrf.QVAPOR( x(i), y(j), z, tWRF ) + wrf.QCLOUD( x(i), y(j), z, tWRF ) + wrf.QRAIN( x(i), y(j), z, tWRF ); % Total water mixing ratio [kg / kg]
				cv.spatial.ql( i, j, z, t ) = wrf.QCLOUD( x(i), y(j), z, tWRF ); % Cloud liquid water mixing ratio [kg / kg]
				cv.spatial.qv( i, j, z, t ) = wrf.QVAPOR( x(i), y(j), z, tWRF ); % Water vapor mixing ratio [kg / kg]
				cv.spatial.qrain( i, j, z, t ) = wrf.QRAIN( x(i), y(j), z, tWRF ); % Rain mixing ratio [kg / kg]
				cv.spatial.T( i, j, z, t ) = wrf.T_PHY( x(i), y(j), z, tWRF ); % Temperature (T) [K]
				cv.spatial.theta( i, j, z, t ) = wrf.T( x(i), y(j), z, tWRF ) + 300; % Potential temperature: perturbation Theta + T0 (300 K) [K]
				cv.spatial.theta_l( i, j, z, t ) = cv.spatial.theta( i, j, z, t) - Lv/cpAir .* (cv.spatial.theta( i, j, z, t) ./ cv.spatial.T( i, j, z, t)) .* cv.spatial.ql( i, j, z, t); % Liquid water potential temperature [K]
				cv.spatial.exner( i, j, z, t ) = (wrf.PSFC( x(i), y(j), tWRF ) ./ wrf.PB( x(i), y(j), z, tWRF )) .^ (R ./ cpAir); % Exner function [-]
				cv.spatial.u( i, j, z, t ) = 0.5 .* ( wrf.U( x(i), y(j), z, tWRF ) + wrf.U( x(i)+1, y(j), z, tWRF ) ); % x-wind component [m / s]
				cv.spatial.v( i, j, z, t ) = 0.5 .* ( wrf.V( x(i), y(j), z, tWRF ) + wrf.V( x(i), y(j)+1, z, tWRF ) ); % y-wind component [m / s]
				cv.spatial.w( i, j, z, t ) = 0.5 .* ( wrf.W( x(i), y(j), z, tWRF ) + wrf.W( x(i), y(j), z+1, tWRF ) ); % z-wind component [m / s]
				% 				cv.spatial.div( i, j, z, t ) = wrf.DIV( x(i), y(j), z, tWRF ); % Divergence [1 / s]
				if isfield(wrf,'QKE')
					cv.spatial.qke( i, j, z, t ) = wrf.QKE( x(i), y(j), z, tWRF ); % Twice TKE from MYNN [m^2 / s^s]
				end
% 				cv.spatial.cf( i, j, z, t) = wrf.CLDFRA( x(i), y(j), z, tWRF ); % Cloud fraction [-]
				cv.spatial.rho( i, j, z, t ) = ( wrf.PB( x(i), y(j), z, tWRF) + wrf.P( x(i), y(j), z, tWRF) ) ./ ( R .* cv.spatial.T( i, j, z, t ) ); % Density [kg / m^3]
				cv.spatial.pmb( i, j, z, t ) = ( wrf.PB( x(i), y(j), z, tWRF) + wrf.P( x(i), y(j), z, tWRF) ) .* 0.01; % Pressure in millibars
				cv.spatial.qs( i, j, z, t ) = qsat( cv.spatial.T( i, j, z, t), cv.spatial.pressure( i, j, z, t ) );
            end
            cv.spatial.pressurew( i, j, :, t) = interp1( squeeze(cv.znu(:, t)), squeeze(cv.spatial.pressure( i, j, :, t)), squeeze(cv.znw( :, t)));
			cv.spatial.lwp(i, j, t) = trapz( squeeze(wrf.height( x(i), y(j), :, tWRF )), squeeze(cv.spatial.ql( i, j, :, t )) ) .* densityAir; % Liquid water path [kg / m^2]

			% Hybrid method to determine zi, following Talbot et al. (2012)
			% Zi = z[ max( | dT/dz * 10^2 | + | dq/dz * 10^5 | ) ]; with dT/dz > 0, dq/dz < 0
			% 1) Create dT/dz and dq/dz vectors
			dTdz = zeros(1, size(cv.height,3)); dqdz = zeros(1, size(cv.height,3));
			for z = 1:findMatch(squeeze(cv.height(i, j, :, t)), 3000) % Constrain to under 3 km
				
				% Vertical temperature gradient
				if z == 1
					dz = cv.height( i, j, z+1, t ) - cv.height( i, j, z, t ); % z(2) - z(1)
					dTdz(z) = ( 1 ./ dz ) .* ( cv.spatial.T( i, j, z+1, t ) - cv.spatial.T( i, j, z, t ) );
				elseif z == size(cv.height, 3)
					dz = cv.height( i, j, z, t ) - cv.height( i, j, z-1, t ); % z(end) - z(end-1)
					dTdz(z) = ( 1 ./ dz ) .* ( cv.spatial.T( i, j, z, t ) - cv.spatial.T( i, j, z-1, t ) );
				else
					dz = cv.height( i, j, z+1, t ) - cv.height( i, j, z, t ); % z(2) - z(1)
					dTdz(z) = ( 1 ./ dz ) .* ( cv.spatial.T( i, j, z+1, t ) - cv.spatial.T( i, j, z, t ) );
				end
				
				% Vertical moisture gradient
				if z == 1
					dz = cv.height( i, j, z+1, t ) - cv.height( i, j, z, t ); % z(2) - z(1)
					dqdz(z) = ( 1 ./ dz ) .* ( cv.spatial.qt( i, j, z+1, t ) - cv.spatial.qt( i, j, z, t ) );
				elseif z == size(cv.height, 3)
					dz = cv.height( i, j, z, t ) - cv.height( i, j, z-1, t ); % z(end) - z(end-1)
					dqdz(z) = ( 1 ./ dz ) .* ( cv.spatial.qt( i, j, z, t ) - cv.spatial.qt( i, j, z-1, t ) );
				else
					dz = cv.height( i, j, z+1, t ) - cv.height( i, j, z, t ); % z(2) - z(1)
					dqdz(z) = ( 1 ./ dz ) .* ( cv.spatial.qt( i, j, z+1, t ) - cv.spatial.qt( i, j, z, t ) );
				end
				
			end
			%%% End zi detection
			
			% 2) Filter out gradients of incorrect sign
			keep = (dTdz > 0) & (dqdz < 0);
			dTdz( ~keep ) = 0; dqdz( ~keep ) = 0;
			
			% 3) Determine z level of max gradients
			sumGrads = abs(dTdz .* 10^2) + abs(dqdz .* 10^5);
			cv.spatial.zi2(i, j, t) = cv.heightAGL(i, j, findMatch(sumGrads, max(sumGrads)), t);
			
			% XZ method
			[~, ~, ~, cv.spatial.zi(i, j, t), ~, ~] = TMP_Inversion_Strength_Cal_V1(squeeze(cv.spatial.T(i, j, :, t)), squeeze(cv.height(i, j, :, t))./1000, cv.spatial.hgt(i, j, t));
			cv.spatial.zi(i, j, t) = 1000.*cv.spatial.zi(i, j, t) - cv.spatial.hgt(i, j, t); % AGL
			
			%%%%%%%%% CLOUD TOP BASED ZI %%%%%%%%
			if ~isempty(find(squeeze(cv.spatial.ql(i, j, :, t)), 1, 'last'))
				cv.spatial.ctop(i, j, t) = find(squeeze(cv.spatial.ql(i, j, :, t)), 1, 'last');
				cv.spatial.idxInversionBase(i, j, t) = cv.spatial.ctop(i, j, t);
				cv.spatial.zi(i, j, t) = cv.heightAGL(i, j, cv.spatial.idxInversionBase(i, j, t), t);
			else
				cv.spatial.ctop(i, j, t) = 0;
			end
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			% ALWAYS DETERMINED FROM ZI. CHANGE NAMING ABOVE TO CHANGE IDX DETERMINATION. (zi [XZ] <--> zi2 [Talbot hybrid])
			cv.spatial.idxInversionBase(i, j, t) = findMatch( squeeze(cv.heightAGL( i, j, :, t)), cv.spatial.zi(i, j, t) );
			
            % ADD IN Grenier Bretherton (2000) [GB2000] SCHEME
            [pi, mu, cv.spatial.dthl(i, j, t), ~, ~] = findInv(squeeze(cv.spatial.theta_l(i, j, :, t)), squeeze(cv.spatial.pressure(i, j, :, t)), squeeze(cv.spatial.pressurew(i, j, :, t)), squeeze(cv.spatial.idxInversionBase(i, j, t)));
            [cv.spatial.dqt(i, j, t), ~, ~] = findJump(squeeze(cv.spatial.qt(i, j, :, t)), squeeze(cv.spatial.pressure(i, j, :, t)), squeeze(cv.spatial.pressurew(i, j, :, t)), squeeze(cv.spatial.idxInversionBase(i, j, t)), mu, pi);
            
			% Integrate TKE to cloud top
			if isfield(wrf,'QKE')
				cv.spatial.tkeInt(i, j, t) = trapz( squeeze(wrf.height( x(i), y(j), 1:cv.spatial.idxInversionBase(i, j, t), tWRF )), squeeze(cv.spatial.qke( i, j, 1:cv.spatial.idxInversionBase(i, j, t), t )) ); % Integrated TKE [ m^3/s^2 ]
			end
			
			% Cloud base
			firstLiquidGridPoint = find( squeeze(cv.spatial.ql( i, j, :, t)) > 0, 1, 'first' ); % Find first grid point in which QCLOUD > 0
			if ~isempty(firstLiquidGridPoint)
				cv.spatial.idxCloudBase(i, j, t) = firstLiquidGridPoint; % Found QCLOUD
			else
				cv.spatial.idxCloudBase(i, j, t) = cv.spatial.idxInversionBase(i, j, t);
			end
			if ~isnan(cv.spatial.idxCloudBase(i, j, t))
				cv.spatial.zb(i, j, t) = cv.heightAGL(i, j, cv.spatial.idxCloudBase(i, j, t), t);
			else
				cv.spatial.zb(i, j, t) = cv.spatial.zi(i, j, t);
			end
			
			% Streamer RTM net flux
			if ~isempty(streamer)
				topIdx = findMatch( streamer.radHeight.*1000, cv.spatial.zi( i, j, t ) ); % botIdx = findMatch( streamer.radHeight.*1000, cv.spatial.zb( i, j, t ) );
				cv.spatial.streamerNet( i, j, t ) = streamer.radMat( topIdx, t ) - streamer.radMat( nlev, t );
			end
			% WRF computed tendencies
			topIdx = cv.spatial.idxInversionBase(i, j, t);
			totalMass = wrf.MU( x(i), y(j), tWRF ) + wrf.MUB( x(i), y(j), tWRF );
			
			if ~isempty(wrf.RTHRATEN)
% 				decoupledRadTen = wrf.RTHRATEN( x(i), y(j), 1:topIdx, tWRF ) ./ totalMass; % Decouple radiative tendency by dividing by mass [Pa K / s / Pa]
				decoupledRadTen = wrf.RTHRATEN( x(i), y(j), 1:topIdx, tWRF ); % As of 3.9, no longer coupled with mass!
				cv.spatial.radtenwrf( i, j, t ) = sum( squeeze( decoupledRadTen ) .* squeeze( cv.DNW( 1:topIdx, t ) ) ) ./ sum( cv.DNW( 1:topIdx, t ) ); % [K / s]
			else
				cv.spatial.radtenwrf( i, j, t ) = 0;
				warning('RTHRATEN DOES NOT EXIST IN INPUT FILE')
			end
			
			% BL averaged relative humidity
			blRH = squeeze(cv.spatial.qv(i, j, 1:topIdx, t)) ./ squeeze(cv.spatial.qs( i, j, 1:topIdx, t));
			cv.spatial.blRH( i, j, t) = sum( squeeze( blRH ) .* squeeze( cv.DNW( 1:topIdx, t ) ) ) ./ sum( cv.DNW( 1:topIdx, t ) ); % [K / s]
			
			% 			decoupledPBLTen = wrf.RTHBLTEN( x(i), y(j), 1:topIdx, tWRF ) ./ totalMass; % Decouple PBL tendency by dividing by mass [Pa K / s / Pa]
			% 			cv.spatial.pbltenwrf( i, j, t ) = sum( squeeze( decoupledPBLTen ) .* squeeze( wrf.DNW( 1:topIdx, tWRF ) ) ) ./ sum( wrf.DNW( 1:topIdx, tWRF ) ); % [K / s]
			
			% 			decoupledSHTen = wrf.RTHSHTEN( x(i), y(j), 1:topIdx, tWRF ) ./ totalMass; % Decouple shallow cumulus tendency by dividing by mass [Pa K / s / Pa]
			% 			cv.spatial.shcutenwrf( i, j, t ) = sum( squeeze( decoupledSHTen ) .* squeeze( wrf.DNW( 1:topIdx, tWRF ) ) ) ./ sum( wrf.DNW( 1:topIdx, tWRF ) ); % [K / s]
			% 			decoupledCUTen = wrf.RTHCUTEN( x(i), y(j), 1:topIdx, tWRF ) ./ totalMass; % Decouple cumulus tendency by dividing by mass [Pa K / s / Pa]
			% 			cv.spatial.cutenwrf( i, j, t ) = sum( squeeze( decoupledCUTen ) .* squeeze( wrf.DNW( 1:topIdx, tWRF ) ) ) ./ sum( wrf.DNW( 1:topIdx, tWRF ) );
			
			% 			cv.spatial.mptenwrf( i, j, t ) = sum( squeeze( wrf.H_DIABATIC( x(i), y(j), 1:topIdx, tWRF ) ) .* squeeze( wrf.DNW( 1:topIdx, tWRF ) ) ) ./ sum( wrf.DNW( 1:topIdx, tWRF ) ); % MP latent heating already decoupled [K / s]
			
			% PBL Averages
		end
	end
end

%% Compute budget equation terms for each grid point in CV

dt = (wrf.time(2) - wrf.time(1)) .* 24 .* 60 .* 60; % dt in seconds

% Initialize computed storage variables with NaN to avoid 0 values at borders
% cv.spatial.dthl = nan(size(cv.spatial.shf));
% cv.spatial.dqt = nan(size(cv.spatial.shf));
cv.spatial.meanThetaL = nan(size(cv.spatial.shf));
cv.spatial.meanQT = nan(size(cv.spatial.shf));

for t = 1:length(cv.time)
	for i = 2:(size(cv.spatial.u, 1) - 1)
		for j = 2:(size(cv.spatial.u, 2) - 1)
			
			%% Some terms to facilitate computation
			theta_l_zi_minus0 = cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t ) , t ); % Theta_l just below inversion [K] [approx to be at inv]
			qt_zi_minus0 = cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t ) , t ); % Total water mixing just point below inversion [K] [approx to be at inv]
			if cv.spatial.idxInversionBase( i, j, t ) > 1
				theta_l_zi_minus = cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t ) - 1, t ); % Theta_l one grid point below inversion [K]
				qt_zi_minus = cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t ) - 1, t ); % Total water mixing ratio one grid point below inversion [K]
			else
				theta_l_zi_minus = theta_l_zi_minus0;
				qt_zi_minus = qt_zi_minus0;
			end
			
			if ~isempty(streamer), deltaF_rad = -cv.spatial.streamerNet( i, j, t ); end % deltaF_rad from streamer
			
			% Inversion height tendency dzi/dt
			if t == 1
				dzi_dt = ( 1 ./ dt ) .* ( cv.spatial.zi( i, j, t+1 ) - cv.spatial.zi( i, j, t ) );
			elseif t == length(cv.time)
				dzi_dt = ( 1 ./ dt ) .* ( cv.spatial.zi( i, j, t ) - cv.spatial.zi( i, j, t-1 ) );
			else
				dzi_dt = ( 1 ./ (2.*dt) ) .* ( cv.spatial.zi( i, j, t+1 ) - cv.spatial.zi( i, j, t-1) );
			end
			
			% Compute on column basis
			% Use mass-weighted boundary layer averages
			
			% Weighted average
			meanU = sum( squeeze(cv.spatial.u( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i, j, t ), t ) );
			meanV = sum( squeeze(cv.spatial.v( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i, j, t ), t ) );
			meanThetaL = sum( squeeze(cv.spatial.theta_l( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i, j, t ), t ) );
			
			meanThetaLxplus = sum( squeeze(cv.spatial.theta_l( i+1, j, 1:cv.spatial.idxInversionBase( i+1, j, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i+1, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i+1, j, t ), t ) );
			meanThetaLxminus = sum( squeeze(cv.spatial.theta_l( i-1, j, 1:cv.spatial.idxInversionBase( i-1, j, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i-1, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i-1, j, t ), t ) );
			meanThetaLyplus = sum( squeeze(cv.spatial.theta_l( i, j+1, 1:cv.spatial.idxInversionBase( i, j+1, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i, j+1, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i, j+1, t ), t ) );
			meanThetaLyminus = sum( squeeze(cv.spatial.theta_l( i, j-1, 1:cv.spatial.idxInversionBase( i, j-1, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i, j-1, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i, j-1, t ), t ) );
			
			dzi_dx = ( 1 ./ (2.*h) ) .* ( cv.spatial.zi( i+1, j, t) - cv.spatial.zi( i-1, j, t) );
			dzi_dy = ( 1 ./ (2.*h) ) .* ( cv.spatial.zi( i, j+1, t) - cv.spatial.zi( i, j-1, t) );
			
			vDotGradzi = meanU .* dzi_dx + meanV .* dzi_dy; % Old method
			
			% d(meanThetaL)/dx and dy
			dthetal_dx = ( 1 ./ (2.*h) ) .* ( meanThetaLxplus - meanThetaLxminus );
			dthetal_dy = ( 1 ./ (2.*h) ) .* ( meanThetaLyplus - meanThetaLyminus );
			
			% v dot grad(meanThetaL)
			vDotGradThetaL = meanU .* dthetal_dx + meanV .* dthetal_dy;
			
			% Subsidence model (Bellon and Stevens, 2011)
			% Assume zw = zi (length scale of large scale dynamics = PBL depth)
			w0 = D .* cv.spatial.zi( i, j, t ); 
			wsub = w0 .* (1 - exp(-1)); % Subsidence velocity at inversion is ...e^(-z/zw), but z = zw = zi --> e^(-1) [Sign convention is opposite in this case because w_e defined positive]
			
% 			w_e = dzi_dt + vDotGradzi + cv.spatial.w( i, j, cv.spatial.idxInversionBase( i, j, t ), t ); % Mass-based entrainment velocity [m / s]
			w_e = dzi_dt + vDotGradzi + wsub;

			if (t > 1) && (abs( (w_e - cv.budgetHeat.w_e(i, j, t - 1)) ./ cv.budgetHeat.w_e(i, j, t - 1)) > 1)
% 				w_e = cv.budgetHeat.w_e(i, j, t - 1);
				warning('w_e change exceeds 100%')
			end
			
			% Store some stuff for answer-checking
			cv.budgetHeat.w_e(i, j, t) = w_e;
			cv.budgetHeat.dzi_dt(i, j, t) = dzi_dt;
			cv.budgetHeat.vDotGradzi(i, j, t) = vDotGradzi;
% 			cv.budgetHeat.wSub(i, j, t) = cv.spatial.w( i, j, cv.spatial.idxInversionBase( i, j, t ), t );
			cv.budgetHeat.wSub(i, j, t) = wsub;
			cv.budgetHeat.wSubAbove(i, j, t) = cv.spatial.w( i, j, cv.spatial.idxInversionBase( i, j, t ) + 2, t ); % Try to take w from 2 cells above inversion
			
			%% Heat budget solved for dTheta_L/dt
			% Advection of zi and ML-deviation term
			cv.budgetHeat.term1(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* ( vDotGradzi ) .* ( theta_l_zi_minus - meanThetaL ); % 1/zi * (v dot grad zi (theta_l_zi- - theta_l))
			% Advection of theta_l
			cv.budgetHeat.term2(i, j, t) = - vDotGradThetaL; % - (v dot grad theta_l)
			% Entrainment warming
% 			cv.budgetHeat.term3(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* w_e .* ( cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )-1, t) ); % 1/zi * w_e * deltaTheta_l_inv
			cv.budgetHeat.term3(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* w_e .* cv.spatial.dthl( i, j, t); % 1/zi * w_e * deltaTheta_l_inv
            % Surface heating
			cv.budgetHeat.term4(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* cv.spatial.shf( i, j, t ) ./ densityAir ./ cpAir; % 1/zi * w'theta_l'_bar
			% Radiative flux divergence
			if ~isempty(streamer)
				cv.budgetHeat.term5(i, j, t) = - ( 1 ./ (cpAir .* mean( cv.spatial.rho( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t ), 3 ) .* cv.spatial.zi( i, j, t )) ) .* deltaF_rad; % - 1/(cp*rho*zi) * deltaF_rad, where deltaF_rad = F_rad(zi) - F_rad(0)
			else
				cv.budgetHeat.term5(i, j, t) = cv.spatial.radtenwrf( i, j, t ); % Use WRF radiation tendency
			end
			% Rain leaving PBL
			if t == 1
				cv.budgetHeat.term6(i, j, t) = ( Lv ./ (cpAir .* cv.spatial.zi( i, j, t ) ) ) .* (cv.spatial.RainSFC( i, j, t+1 ) - cv.spatial.RainSFC( i, j, t) ) ./ dt ./ 1000 .* densityWater ./ densityAir; % L/(cp*zi) * deltaF_p
			else
				cv.budgetHeat.term6(i, j, t) = ( Lv ./ (cpAir .* cv.spatial.zi( i, j, t ) ) ) .* (cv.spatial.RainSFC( i, j, t ) - cv.spatial.RainSFC( i, j, t-1 )) ./ dt ./ 1000 .* densityWater ./ densityAir; % L/(cp*zi) * deltaF_p
			end
			
			%% Moisture budget solved for dqt/dt
% 			meanqt = mean( cv.spatial.qt( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t), 3 );
			meanqt = sum( squeeze(cv.spatial.qt( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i, j, t ), t ) );
			
			meanqtxplus = sum( squeeze(cv.spatial.qt( i+1, j, 1:cv.spatial.idxInversionBase( i+1, j, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i+1, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i+1, j, t ), t ) );
			meanqtxminus = sum( squeeze(cv.spatial.qt( i-1, j, 1:cv.spatial.idxInversionBase( i-1, j, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i-1, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i-1, j, t ), t ) );
			meanqtyplus = sum( squeeze(cv.spatial.qt( i, j+1, 1:cv.spatial.idxInversionBase( i, j+1, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i, j+1, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i, j+1, t ), t ) );
			meanqtyminus = sum( squeeze(cv.spatial.qt( i, j-1, 1:cv.spatial.idxInversionBase( i, j-1, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i, j-1, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i, j-1, t ), t ) );
			
			dqt_dx = ( 1 ./ (2.*h) ) .* ( meanqtxplus - meanqtxminus );
			dqt_dy = ( 1 ./ (2.*h) ) .* ( meanqtyplus - meanqtyminus );

			vDotGradqt = meanU .* dqt_dx + meanV .* dqt_dy;
			
			% Collect for advection sounding if requested
			if makeSounding
				thetaLAdvection(t) = vDotGradThetaL;
				qtAdvection(t) = vDotGradqt;
				ziAdvection(t) = vDotGradzi;
				wsub(t) = cv.budgetHeat.wSub(i, j, t);
				wsubAbove(t) = cv.budgetHeat.wSubAbove(i, j, t);
			end
			
			% TEMP: SAVE ADVECTION TERMS SEPARATELY
			cv.temp.meanU(i, j, t) = meanU;
			cv.temp.meanV(i, j, t) = meanV;
			cv.temp.dthetal_dx(i, j, t) = dthetal_dx;
			cv.temp.dthetal_dy(i, j, t) = dthetal_dy;
% 			cv.temp.dthetal_inv(i, j, t) = cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )-1, t);
			cv.temp.dthetal_inv2(i, j, t) = cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - meanThetaL;
			cv.temp.meanThetaL(i, j, t) = meanThetaL;
			cv.temp.theta_l_inv(i, j, t) = cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t);
			cv.temp.dqt_dx(i, j, t) = dqt_dx;
			cv.temp.dqt_dy(i, j, t) = dqt_dy;
% 			cv.temp.dqt_inv(i, j, t) = cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )-1, t);
			cv.temp.dqt_inv2(i, j, t) = cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - meanqt;
			cv.temp.meanqt(i, j, t) = meanqt;
			cv.temp.qt_inv(i, j, t) = cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t);
			
			% Advection of zi and ML-deviation term
			cv.budgetMoisture.term1(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* vDotGradzi .* ( qt_zi_minus - meanqt ); % (1 / zi) .* vDotGradzi .* ( q_t_zi - q_t_BL )
			% Advection of q_t
			cv.budgetMoisture.term2(i, j, t) = - vDotGradqt; % - vDotGradqt
			% Entrainment drying
% 			cv.budgetMoisture.term3(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* w_e .* ( cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )-1, t) ); % w_e / zi .* dqt_inv
			cv.budgetMoisture.term3(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* w_e .* cv.spatial.dqt(i, j, t); % w_e / zi .* dqt_inv
            % Surface moistening
% 			cv.budgetMoisture.term4(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* cv.spatial.lhf( i, j, t ) ./ Lv ./ densityWater; %
			cv.budgetMoisture.term4(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* cv.spatial.lhf( i, j, t ) ./ Lv ./ 1.2; % HY TEMPORARY CHECK WITH DENSITY OF AIR RATHER THAN WATER
			% 			cv.budgetMoisture.term4(i-1, j-1, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* cv.spatial.qfx( i, j, t ) ./ densityWater; % Use QFX rather than LHF (e.g. when using RUC)
			% Rain leaving PBL
			if t == 1
				cv.budgetMoisture.term5(i, j, t) = - ( 1 ./ cv.spatial.zi( i, j, t ) ) .* (cv.spatial.RainSFC( i, j, t+1 ) - cv.spatial.RainSFC( i, j, t) ) ./ dt ./ 1000 .* densityWater ./ densityAir; % rain rate is negative
			else
				cv.budgetMoisture.term5(i, j, t) = - ( 1 ./ cv.spatial.zi( i, j, t ) ) .* (cv.spatial.RainSFC( i, j, t ) - cv.spatial.RainSFC( i, j, t-1 )) ./ dt ./ 1000 .* densityWater ./ densityAir; % rain rate is negative
			end
			% Compute heating rates and moisture tendency using centered finite difference (forward and backward at boundaries)
			if t == 1
				cv.budgetHeat.dTheta_L_dt_FD( i, j, t ) = (1 ./ dt ) .* ( mean( cv.spatial.theta_l( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t+1), 3 ) - mean( cv.spatial.theta_l( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t ), 3 ) );
				cv.budgetHeat.dTheta_dt_FD( i, j, t ) = (1 ./ dt ) .* ( mean( cv.spatial.theta( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t+1), 3 ) - mean( cv.spatial.theta( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t ), 3 ) );
				cv.budgetMoisture.dqt_dt_FD( i, j, t ) = (1 ./ dt ) .* ( mean( cv.spatial.qt( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t+1), 3 ) - mean( cv.spatial.qt( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t ), 3 ) );
			elseif t == length(cv.time)
				cv.budgetHeat.dTheta_L_dt_FD( i, j, t ) = (1 ./ dt ) .* ( mean( cv.spatial.theta_l( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t), 3 ) - mean( cv.spatial.theta_l( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t-1 ), 3 ) );
				cv.budgetHeat.dTheta_dt_FD( i, j, t ) = (1 ./ dt ) .* ( mean( cv.spatial.theta( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t), 3 ) - mean( cv.spatial.theta( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t-1 ), 3 ) );
				cv.budgetMoisture.dqt_dt_FD( i, j, t ) = (1 ./ dt ) .* ( mean( cv.spatial.qt( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t), 3 ) - mean( cv.spatial.qt( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t-1 ), 3 ) );
			else
				cv.budgetHeat.dTheta_L_dt_FD( i, j, t ) = (1 ./ (2.*dt) ) .* ( mean( cv.spatial.theta_l( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t+1), 3 ) - mean( cv.spatial.theta_l( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t-1 ), 3 ) );
				cv.budgetHeat.dTheta_dt_FD( i, j, t ) = (1 ./ (2.*dt) ) .* ( mean( cv.spatial.theta( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t+1), 3 ) - mean( cv.spatial.theta( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t-1 ), 3 ) );
				cv.budgetMoisture.dqt_dt_FD( i, j, t ) = (1 ./ (2.*dt) ) .* ( mean( cv.spatial.qt( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t+1), 3 ) - mean( cv.spatial.qt( i, j, 1:cv.spatial.idxInversionBase( i, j, t ), t-1 ), 3 ) );
			end
			
			% Cloud base tendencies
			Tcb = cv.spatial.T(i, j, cv.spatial.idxCloudBase(i, j, t), t); % Temperature T at cloud base
			cv.cloudBase.dzb_dqt(i, j, t) = (R .* Tcb) ./ (g .* meanqt) .* (1 - (Lv.*R)./(cpAir.*Rv.*Tcb)).^(-1); % Change in cloud base due to change in q_t (i.e. efficiency)
			cv.cloudBase.dzb_dthetal(i, j, t) = 1 ./ ( g ./ (cpAir .* cv.spatial.exner(i, j, cv.spatial.idxCloudBase(i, j, t), t)) .* (1 - (cpAir.*Rv.*Tcb)./(R.*Lv)) ); % Change in cloud base due to change in theta_l (i.e. efficiency)
			% Compute cloud base tendency using centered finite difference (forward and backward at boundaries)
			if t == 1
				cv.cloudBase.dzb_dt_FD(i, j, t) = (1 ./ dt ) .* ( cv.spatial.zb(i, j, t+1) - cv.spatial.zb(i, j, t) );
			elseif t == length(cv.time)
				cv.cloudBase.dzb_dt_FD(i, j, t) = (1 ./ dt ) .* ( cv.spatial.zb(i, j, t) - cv.spatial.zb(i, j, t-1) );
			else
				cv.cloudBase.dzb_dt_FD(i, j, t) = (1 ./ (2.*dt) ) .* ( cv.spatial.zb(i, j, t+1) - cv.spatial.zb(i, j, t-1) );
			end
			
			cv.temp.qtDeviation(i, j, t) = qt_zi_minus - meanqt;
			cv.temp.tlDeviation(i, j, t) = theta_l_zi_minus - meanThetaL;
			cv.temp.vDotGradzi(i, j, t) = vDotGradzi;
			
			% Store calculated variables for output
% 			cv.spatial.dthl(i, j, t) = cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )-1, t);
% 			cv.spatial.dqt(i, j, t) = cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )-1, t);
			cv.spatial.meanThetaL(i, j, t ) = meanThetaL;
			cv.spatial.meanQT(i,j,t) = meanqt;
			
		end
	end
end

%% Sum up to obtain heating rate and moisture tendency!
cv.budgetHeat.dTheta_L_dt = cv.budgetHeat.term1 + cv.budgetHeat.term2 + cv.budgetHeat.term3 + cv.budgetHeat.term4 + cv.budgetHeat.term5 + cv.budgetHeat.term6;
cv.budgetMoisture.dqt_dt = cv.budgetMoisture.term1 + cv.budgetMoisture.term2 + cv.budgetMoisture.term3 + cv.budgetMoisture.term4 + cv.budgetMoisture.term5;

% Make assumption that residual between heating rate/moisture tendency is solely due to the entrainment term (i.e. all other terms are correct)
% --> back out entrainment term
cv.budgetHeat.term3closure = cv.budgetHeat.dTheta_L_dt_FD - cv.budgetHeat.term1 - cv.budgetHeat.term2 - cv.budgetHeat.term4 - cv.budgetHeat.term5 - cv.budgetHeat.term6;
cv.budgetMoisture.term3closure = cv.budgetMoisture.dqt_dt_FD - cv.budgetMoisture.term1 - cv.budgetMoisture.term2 - cv.budgetMoisture.term4 - cv.budgetMoisture.term5;

%% Use FD HR to back out w_e from HB and MB residuals
% Use these to verify residual only contains entrainment term
for t = 1:length(cv.time)
	for i = 2:(size(cv.spatial.u, 1)-1)
		for j = 2:(size(cv.spatial.u, 2)-1)
% 			cv.budgetMoisture.backedOutw_e(i, j, t) = ( cv.budgetMoisture.dqt_dt_FD(i, j, t) - cv.budgetMoisture.term1(i, j, t) - cv.budgetMoisture.term2(i, j, t) - cv.budgetMoisture.term4(i, j, t) - cv.budgetMoisture.term5(i, j, t) ) .* cv.spatial.zi( i, j, t ) ./ ( cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )-1, t) ); % dqt_dt .* z_i ./ dqt_inv
			cv.budgetMoisture.backedOutw_e(i, j, t) = ( cv.budgetMoisture.dqt_dt_FD(i, j, t) - cv.budgetMoisture.term1(i, j, t) - cv.budgetMoisture.term2(i, j, t) - cv.budgetMoisture.term4(i, j, t) - cv.budgetMoisture.term5(i, j, t) ) .* cv.spatial.zi( i, j, t ) ./ ( cv.spatial.dqt(i, j, t) ); % dqt_dt .* z_i ./ dqt_inv
% 			cv.budgetHeat.term3FromMoistureBudget(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* cv.budgetMoisture.backedOutw_e(i, j, t) .* ( cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )-1, t) );
			cv.budgetHeat.term3FromMoistureBudget(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* cv.budgetMoisture.backedOutw_e(i, j, t) .* ( cv.spatial.dthl( i, j, t) );
			
% 			cv.budgetHeat.backedOutw_e(i, j, t) = ( cv.budgetHeat.dTheta_L_dt_FD(i, j, t) - cv.budgetHeat.term1(i, j, t) - cv.budgetHeat.term2(i, j, t) - cv.budgetHeat.term4(i, j, t) - cv.budgetHeat.term5(i, j, t) - cv.budgetHeat.term6(i, j, t) ) .* cv.spatial.zi( i, j, t ) ./ ( cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )-1, t) );
cv.budgetHeat.backedOutw_e(i, j, t) = ( cv.budgetHeat.dTheta_L_dt_FD(i, j, t) - cv.budgetHeat.term1(i, j, t) - cv.budgetHeat.term2(i, j, t) - cv.budgetHeat.term4(i, j, t) - cv.budgetHeat.term5(i, j, t) - cv.budgetHeat.term6(i, j, t) ) .* cv.spatial.zi( i, j, t ) ./ ( cv.spatial.dthl( i, j, t) );
% 			cv.budgetHeat.term3FromHeatBudget(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* cv.budgetHeat.backedOutw_e(i, j, t) .* ( cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - cv.spatial.theta_l( i, j, cv.spatial.idxInversionBase( i, j, t )-1, t) );
			cv.budgetHeat.term3FromHeatBudget(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* cv.budgetHeat.backedOutw_e(i, j, t) .* ( cv.spatial.dthl( i, j, t) );
% 			cv.budgetMoisture.term3FromHeatBudget(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* cv.budgetHeat.backedOutw_e(i, j, t) .* ( cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )-1, t) ); % w_e / zi .* dqt_inv
			cv.budgetMoisture.term3FromHeatBudget(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* cv.budgetHeat.backedOutw_e(i, j, t) .* ( cv.spatial.dqt(i, j, t) ); % w_e / zi .* dqt_inv
			
			cv.budgetHeat.meanw_e(i, j, t) = mean( [cv.budgetHeat.backedOutw_e(i, j, t), cv.budgetMoisture.backedOutw_e(i, j, t)] );
			cv.budgetHeat.term3mean(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* cv.budgetHeat.meanw_e(i, j, t) .* ( cv.spatial.dthl( i, j, t) );
% 			cv.budgetMoisture.term3mean(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* cv.budgetHeat.meanw_e(i, j, t) .* ( cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )+1, t) - cv.spatial.qt( i, j, cv.spatial.idxInversionBase( i, j, t )-1, t) ); % w_e / zi .* dqt_inv
			cv.budgetMoisture.term3mean(i, j, t) = ( 1 ./ cv.spatial.zi( i, j, t ) ) .* cv.budgetHeat.meanw_e(i, j, t) .* ( cv.spatial.dqt(i, j, t) ); % w_e / zi .* dqt_inv
		end
	end
end

% Sum up
cv.budgetMoisture.dqt_dt_HB = cv.budgetMoisture.term1 + cv.budgetMoisture.term2 + cv.budgetMoisture.term3FromHeatBudget + cv.budgetMoisture.term4 + cv.budgetMoisture.term5;
cv.budgetMoisture.dqt_dt_mean = cv.budgetMoisture.term1 + cv.budgetMoisture.term2 + cv.budgetMoisture.term3mean + cv.budgetMoisture.term4 + cv.budgetMoisture.term5;
cv.budgetHeat.dTheta_L_dt_MB = cv.budgetHeat.term1 + cv.budgetHeat.term2 + cv.budgetHeat.term3FromMoistureBudget + cv.budgetHeat.term4 + cv.budgetHeat.term5 + cv.budgetHeat.term6;
cv.budgetHeat.dTheta_L_dt_HB = cv.budgetHeat.term1 + cv.budgetHeat.term2 + cv.budgetHeat.term3FromHeatBudget + cv.budgetHeat.term4 + cv.budgetHeat.term5 + cv.budgetHeat.term6;
cv.budgetHeat.dTheta_L_dt_mean = cv.budgetHeat.term1 + cv.budgetHeat.term2 + cv.budgetHeat.term3mean + cv.budgetHeat.term4 + cv.budgetHeat.term5 + cv.budgetHeat.term6;

%% Cloud base tendency
% dzb/dt = dzb/dthetal .* dthetal/dt + dzb/dqt .* dqt/dt;
cv.cloudBase.dzb_dt_heat = cv.cloudBase.dzb_dthetal .* cv.budgetHeat.dTheta_L_dt_FD;
cv.cloudBase.dzb_dt_moisture = cv.cloudBase.dzb_dqt .* cv.budgetMoisture.dqt_dt_FD;
cv.cloudBase.dzb_dt = cv.cloudBase.dzb_dt_heat + cv.cloudBase.dzb_dt_moisture; % Use FD heating rate/moisture tendency

% Separate
cv.cloudBase.dzb_dt_h1 = cv.cloudBase.dzb_dthetal .* cv.budgetHeat.term1;
cv.cloudBase.dzb_dt_h2 = cv.cloudBase.dzb_dthetal .* cv.budgetHeat.term2;
cv.cloudBase.dzb_dt_h3 = cv.cloudBase.dzb_dthetal .* cv.budgetHeat.term3mean;
cv.cloudBase.dzb_dt_h4 = cv.cloudBase.dzb_dthetal .* cv.budgetHeat.term4;
cv.cloudBase.dzb_dt_h5 = cv.cloudBase.dzb_dthetal .* cv.budgetHeat.term5;
cv.cloudBase.dzb_dt_h6 = cv.cloudBase.dzb_dthetal .* cv.budgetHeat.term6;

cv.cloudBase.dzb_dt_m1 = cv.cloudBase.dzb_dqt .* cv.budgetMoisture.term1;
cv.cloudBase.dzb_dt_m2 = cv.cloudBase.dzb_dqt .* cv.budgetMoisture.term2;
cv.cloudBase.dzb_dt_m3 = cv.cloudBase.dzb_dqt .* cv.budgetMoisture.term3mean;
cv.cloudBase.dzb_dt_m4 = cv.cloudBase.dzb_dqt .* cv.budgetMoisture.term4;
cv.cloudBase.dzb_dt_m5 = cv.cloudBase.dzb_dqt .* cv.budgetMoisture.term5;

for i = 2:(size(cv.spatial.u, 1)-1)
	for j = 2:(size(cv.spatial.u, 2)-1)
		cv.cloudBase.wrfThickness(i, j, 1) = cv.spatial.zi(i, j, 1) - cv.spatial.zb(i, j, 1);
		cv.cloudBase.predictedThickness(i, j, 1) = cv.cloudBase.wrfThickness(i, j, 1);
		cv.cloudBase.predictedzb(i, j, 1) = cv.spatial.zb(i, j, 2);
		cv.cloudBase.predictedzi(i, j, 1) = cv.spatial.zi(i, j, 2);
		for t = 2:length(cv.time)
			cv.cloudBase.predictedzi(i, j, t) = cv.cloudBase.predictedzi(i, j, t-1) + (cv.budgetHeat.meanw_e(i, j, t-1) - cv.budgetHeat.wSub(i, j, t-1)) .* dt; % dzi/dt + adv = w_e + w_sub
			cv.cloudBase.predictedzb(i, j, t) = cv.cloudBase.predictedzb(i, j, t-1) + cv.cloudBase.dzb_dt(i, j, t-1) .* dt;
			if cv.cloudBase.predictedzb(i, j, t) < 0, cv.cloudBase.predictedzb(i, j, t) = 0; end
			cv.cloudBase.predictedThickness(i, j, t) = cv.spatial.zi(i, j, t) - cv.cloudBase.predictedzb(i, j, t);
			if cv.cloudBase.predictedThickness(i, j, t) < 0, cv.cloudBase.predictedThickness(i, j, t) = 0; end
			cv.cloudBase.wrfThickness(i, j, t) = cv.spatial.zi(i, j, t) - cv.spatial.zb(i, j, t);
			if cv.cloudBase.wrfThickness(i, j, t) < 0, cv.cloudBase.wrfThickness(i, j, t) = 0; end
		end
		
	end
end

%% Spatially average select variables. Should you want more, just add them here into the correct cell array.
cv.var3D = {'shf', 'lhf', 'psfc', 'hgt', 'PBLH', 'lwp', 'dthl', 'dqt', 'meanThetaL', 'meanQT', 'idxInversionBase', 'idxCloudBase', 'zi', 'zb', 'RainSFC', 'radtenwrf'};
cv.var4D = {'qt', 'theta', 'rho', 'u', 'v', 'w', 'T', 'ql', 'theta_l', 'exner'}; % ,'cf'

for var = 1:length(cv.var3D)
	try
	% Take average of everything except boundaries
	temp = cv.spatial.(cv.var3D{var})(2:end-1, 2:end-1, :);
	cv.avg.(cv.var3D{var}) = nanmean( nanmean(temp, 1), 2);
	cv.avg.(cv.var3D{var}) = squeeze(cv.avg.(cv.var3D{var}));
	catch
	end
end

for var = 1:length(cv.var4D)
	% Take average of everything except boundaries
	try
	temp = cv.spatial.(cv.var4D{var})(2:end-1, 2:end-1, :, :);
	cv.avg.(cv.var4D{var})(:, :) = nanmean( nanmean( temp, 1), 2);
	catch
	end
end

%% DO EVERYTHING BUT ON CV AVERAGE RATHER THAN COLUMN BY COLUMN

cv.avg.idxInversionBase = round(cv.avg.idxInversionBase);
cv.avg.idxCloudBase = round(cv.avg.idxCloudBase);
cv.avg.zb = round(cv.avg.zb);

for t = 1:length(cv.time)
	
	%% Some terms to facilitate computation
	theta_l_zi_minus0 = cv.avg.theta_l( cv.avg.idxInversionBase( t ) , t ); % Theta_l just below inversion [K] [approx to be at inv]
	qt_zi_minus0 = cv.avg.qt( cv.avg.idxInversionBase( t ) , t ); % Total water mixing just point below inversion [K] [approx to be at inv]
	if cv.avg.idxInversionBase( t ) < 1
		theta_l_zi_minus = cv.avg.theta_l( cv.avg.idxInversionBase( t ) - 1, t ); % Theta_l one grid point below inversion [K]
		qt_zi_minus = cv.avg.qt( cv.avg.idxInversionBase( t ) - 1, t ); % Total water mixing ratio one grid point below inversion [K]
	else
		theta_l_zi_minus = theta_l_zi_minus0;
		qt_zi_minus = qt_zi_minus0;
	end
	if ~isempty(streamer), deltaF_rad = -cv.avg.streamerNet( t ); end % deltaF_rad from streamer
	
	if t == 1
		dzi_dt = ( 1 ./ dt ) .* ( cv.avg.zi( t+1 ) - cv.avg.zi( t ) );
	elseif t == length(cv.time)
		dzi_dt = ( 1 ./ dt ) .* ( cv.avg.zi( t ) - cv.avg.zi( t-1 ) );
	else
		dzi_dt = ( 1 ./ (2.*dt) ) .* ( cv.avg.zi( t+1 ) - cv.avg.zi( t-1) );
	end
	
	% Compute on column basis
	% Just use mean values for the moment..?
	
	% Weighted average
	meanU = sum( squeeze(cv.avg.u( 1:cv.avg.idxInversionBase( t ), t)) .* squeeze( cv.DNW( 1:cv.avg.idxInversionBase( t ), t ) ) ) ./ sum( cv.DNW( 1:cv.avg.idxInversionBase( t ), t ) );
	meanV = sum( squeeze(cv.avg.v( 1:cv.avg.idxInversionBase( t ), t)) .* squeeze( cv.DNW( 1:cv.avg.idxInversionBase(  t ), t ) ) ) ./ sum( cv.DNW( 1:cv.avg.idxInversionBase(  t ), t ) );
	meanThetaL = sum( squeeze(cv.avg.theta_l(  1:cv.avg.idxInversionBase(  t ), t)) .* squeeze( cv.DNW( 1:cv.avg.idxInversionBase(  t ), t ) ) ) ./ sum( cv.DNW( 1:cv.avg.idxInversionBase(  t ), t ) );
	
	% U(x) direction
	meanThetaLxminus = 0;
	meanThetaLxplus = 0;
	for j = 1:size(cv.spatial.u, 2) % Yes, the indices are correct [they go normal to the velocity]
		meanThetaLxminus = meanThetaLxminus + sum( squeeze(cv.spatial.theta_l( 1, j, 1:cv.spatial.idxInversionBase( 1, j, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( 1, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( 1, j, t ), t ) );
		meanThetaLxplus = meanThetaLxplus + sum( squeeze(cv.spatial.theta_l( end, j, 1:cv.spatial.idxInversionBase( end, j, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( end, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( end, j, t ), t ) );
	end
	meanThetaLxminus = meanThetaLxminus ./ size(cv.spatial.u, 2);
	meanThetaLxplus = meanThetaLxplus ./ size(cv.spatial.u, 2);
	% V(y) direction
	meanThetaLyminus = 0;
	meanThetaLyplus = 0;
	for i = 1:size(cv.spatial.v, 1)
		meanThetaLyminus = meanThetaLyminus + sum( squeeze(cv.spatial.theta_l( i, 1, 1:cv.spatial.idxInversionBase( i, 1, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i, 1, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i, 1, t ), t ) );
		meanThetaLyplus = meanThetaLyplus + sum( squeeze(cv.spatial.theta_l( i, end, 1:cv.spatial.idxInversionBase( i, end, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i, end, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i, end, t ), t ) );
	end
	meanThetaLyminus = meanThetaLyminus ./ size(cv.spatial.v, 1);
	meanThetaLyplus = meanThetaLyplus ./ size(cv.spatial.v, 1);
	% 			meanThetaLxplus = sum( squeeze(cv.avg.theta_l( i+1, j, 1:cv.avg.idxInversionBase( i+1, j, t ), t)) .* squeeze( cv.DNW( 1:cv.avg.idxInversionBase( i+1, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.avg.idxInversionBase( i+1, j, t ), t ) );
	%			meanThetaLxminus = sum( squeeze(cv.avg.theta_l( i-1, j, 1:cv.avg.idxInversionBase( i-1, j, t ), t)) .* squeeze( cv.DNW( 1:cv.avg.idxInversionBase( i-1, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.avg.idxInversionBase( i-1, j, t ), t ) );
	%			meanThetaLyplus = sum( squeeze(cv.avg.theta_l( i, j+1, 1:cv.avg.idxInversionBase( i, j+1, t ), t)) .* squeeze( cv.DNW( 1:cv.avg.idxInversionBase( i, j+1, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.avg.idxInversionBase( i, j+1, t ), t ) );
	%			meanThetaLyminus = sum( squeeze(cv.avg.theta_l( i, j-1, 1:cv.avg.idxInversionBase( i, j-1, t ), t)) .* squeeze( cv.DNW( 1:cv.avg.idxInversionBase( i, j-1, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.avg.idxInversionBase( i, j-1, t ), t ) );
	
	%			dzi_dx = ( 1 ./ (2.*h) ) .* ( cv.avg.zi( i+1, j, t) - cv.avg.zi( i-1, j, t) );
	%			dzi_dy = ( 1 ./ (2.*h) ) .* ( cv.avg.zi( i, j+1, t) - cv.avg.zi( i, j-1, t) );
	
	%			vDotGradzi = meanU .* dzi_dx + meanV .* dzi_dy; % Old method
	vDotGradzi = 0; % This term is typically 0 unless there is a step change from column-to-column (rare), but causes a lot of instability. Zeroing out this term for stability.
	
	% d(meanThetaL)/dx and dy
	
	dthetal_dx = ( 1 ./ ((length(x)-1).*h) ) .* ( meanThetaLxplus - meanThetaLxminus );
	dthetal_dy = ( 1 ./ ((length(y)-1).*h) ) .* ( meanThetaLyplus - meanThetaLyminus );
	
	% v dot grad(meanThetaL)
	vDotGradThetaL = meanU .* dthetal_dx + meanV .* dthetal_dy;
	
	% Subsidence model (Bellon and Stevens, 2011)
	% Assume zw = zi (length scale of large scale dynamics = PBL depth)
	w0 = D .* cv.avg.zi(  t );
	wsub = w0 .* (1 - exp(-1)); % Subsidence velocity at inversion is ...e^(-z/zw), but z = zw = zi --> e^(-1) [Sign convention is opposite in this case because w_e defined positive]
	
	% 			w_e = dzi_dt + vDotGradzi + cv.avg.w(  cv.avg.idxInversionBase(  t ), t ); % Mass-based entrainment velocity [m / s]
	w_e = dzi_dt + vDotGradzi + wsub;
	
	if (t > 1) && (abs( (w_e - cv.avg.budgetHeat.w_e( t - 1)) ./ cv.avg.budgetHeat.w_e( t - 1)) > 1)
		% 				w_e = cv.budgetHeat.w_e( t - 1);
		disp('change too big')
	end
	
	% Store some stuff for answer-checking
	cv.avg.budgetHeat.w_e( t) = w_e;
	cv.avg.budgetHeat.dzi_dt( t) = dzi_dt;
	cv.avg.budgetHeat.vDotGradzi( t) = vDotGradzi;
	% 			cv.avg.budgetHeat.wSub( t) = cv.avg.w(  cv.avg.idxInversionBase(  t ), t );
	cv.avg.budgetHeat.wSub( t) = wsub;
	cv.avg.budgetHeat.wSubAbove( t) = cv.avg.w(  cv.avg.idxInversionBase(  t ) + 2, t ); % Try to take w from 2 cells above inversion
	
	%% Heat budget solved for dTheta_L/dt
	% Advection of zi and ML-deviation term
	cv.avg.budgetHeat.term1( t) = ( 1 ./ cv.avg.zi(  t ) ) .* ( vDotGradzi ) .* ( theta_l_zi_minus - meanThetaL ); % 1/zi * (v dot grad zi (theta_l_zi- - theta_l))
	% Advection of theta_l
	cv.avg.budgetHeat.term2( t) = - vDotGradThetaL; % - (v dot grad theta_l)
	% Entrainment warming
% 	cv.avg.budgetHeat.term3( t) = ( 1 ./ cv.avg.zi(  t ) ) .* w_e .* ( cv.avg.theta_l(  cv.avg.idxInversionBase(  t )+1, t) - cv.avg.theta_l(  cv.avg.idxInversionBase(  t )-1, t) ); % 1/zi * w_e * deltaTheta_l_inv
	cv.avg.budgetHeat.term3( t) = ( 1 ./ cv.avg.zi(  t ) ) .* w_e .* ( cv.avg.dthl( t ) ); % 1/zi * w_e * deltaTheta_l_inv
	% Surface heating
	cv.avg.budgetHeat.term4( t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.shf(  t ) ./ densityAir ./ cpAir; % 1/zi * w'theta_l'_bar
	% Radiative flux divergence
	if ~isempty(streamer)
		cv.avg.budgetHeat.term5( t) = - ( 1 ./ (cpAir .* mean( cv.avg.rho(  1:cv.avg.idxInversionBase(  t ), t ), 3 ) .* cv.avg.zi(  t )) ) .* deltaF_rad; % - 1/(cp*rho*zi) * deltaF_rad, where deltaF_rad = F_rad(zi) - F_rad(0)
	else
		cv.avg.budgetHeat.term5( t) = cv.avg.radtenwrf(  t ); % Use WRF radiation tendency
	end
	% Rain leaving PBL
	if t == 1
		cv.avg.budgetHeat.term6( t) = ( Lv ./ (cpAir .* cv.avg.zi(  t ) ) ) .* (cv.avg.RainSFC(  t+1 ) - cv.avg.RainSFC(  t) ) ./ dt ./ 1000 .* densityWater ./ densityAir; % L/(cp*zi) * deltaF_p
	else
		cv.avg.budgetHeat.term6( t) = ( Lv ./ (cpAir .* cv.avg.zi(  t ) ) ) .* (cv.avg.RainSFC(  t ) - cv.avg.RainSFC(  t-1 )) ./ dt ./ 1000 .* densityWater ./ densityAir; % L/(cp*zi) * deltaF_p
	end
	
	%% Moisture budget solved for dqt/dt
	% 			meanqt = mean( cv.avg.qt(  1:cv.avg.idxInversionBase(  t ), t), 3 );
	meanqt = sum( squeeze(cv.avg.qt(  1:cv.avg.idxInversionBase(  t ), t)) .* squeeze( cv.DNW( 1:cv.avg.idxInversionBase(  t ), t ) ) ) ./ sum( cv.DNW( 1:cv.avg.idxInversionBase(  t ), t ) );
	
	% U(x) direction
	meanqtxminus = 0;
	meanqtxplus = 0;
	for j = 1:size(cv.spatial.u, 2) % Yes, the indices are correct [they go normal to the velocity]
		meanqtxminus = meanqtxminus + sum( squeeze(cv.spatial.qt( 1, j, 1:cv.spatial.idxInversionBase( 1, j, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( 1, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( 1, j, t ), t ) );
		meanqtxplus = meanqtxplus + sum( squeeze(cv.spatial.qt( end, j, 1:cv.spatial.idxInversionBase( end, j, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( end, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( end, j, t ), t ) );
	end
	meanqtxminus = meanqtxminus ./ size(cv.spatial.u, 2);
	meanqtxplus = meanqtxplus ./ size(cv.spatial.u, 2);
	% V(y) direction
	meanqtyminus = 0;
	meanqtyplus = 0;
	for i = 1:size(cv.spatial.v, 1)
		meanqtyminus = meanqtyminus + sum( squeeze(cv.spatial.qt( i, 1, 1:cv.spatial.idxInversionBase( i, 1, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i, 1, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i, 1, t ), t ) );
		meanqtyplus = meanqtyplus + sum( squeeze(cv.spatial.qt( i, end, 1:cv.spatial.idxInversionBase( i, end, t ), t)) .* squeeze( cv.DNW( 1:cv.spatial.idxInversionBase( i, end, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.spatial.idxInversionBase( i, end, t ), t ) );
	end
	meanqtyminus = meanqtyminus ./ size(cv.spatial.v, 1);
	meanqtyplus = meanqtyplus ./ size(cv.spatial.v, 1);
	%			meanqtxplus = sum( squeeze(cv.avg.qt( i+1, j, 1:cv.avg.idxInversionBase( i+1, j, t ), t)) .* squeeze( cv.DNW( 1:cv.avg.idxInversionBase( i+1, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.avg.idxInversionBase( i+1, j, t ), t ) );
	%			meanqtxminus = sum( squeeze(cv.avg.qt( i-1, j, 1:cv.avg.idxInversionBase( i-1, j, t ), t)) .* squeeze( cv.DNW( 1:cv.avg.idxInversionBase( i-1, j, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.avg.idxInversionBase( i-1, j, t ), t ) );
	%			meanqtyplus = sum( squeeze(cv.avg.qt( i, j+1, 1:cv.avg.idxInversionBase( i, j+1, t ), t)) .* squeeze( cv.DNW( 1:cv.avg.idxInversionBase( i, j+1, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.avg.idxInversionBase( i, j+1, t ), t ) );
	%			meanqtyminus = sum( squeeze(cv.avg.qt( i, j-1, 1:cv.avg.idxInversionBase( i, j-1, t ), t)) .* squeeze( cv.DNW( 1:cv.avg.idxInversionBase( i, j-1, t ), t ) ) ) ./ sum( cv.DNW( 1:cv.avg.idxInversionBase( i, j-1, t ), t ) );
	
	dqt_dx = ( 1 ./ ((length(x)-1).*h) ) .* ( meanqtxplus - meanqtxminus );
	dqt_dy = ( 1 ./ ((length(y)-1).*h) ) .* ( meanqtyplus - meanqtyminus );
	
	vDotGradqt = meanU .* dqt_dx + meanV .* dqt_dy;
	
	% Collect for advection sounding if requested
	if makeSounding
		thetaLAdvection(t) = vDotGradThetaL;
		qtAdvection(t) = vDotGradqt;
		ziAdvection(t) = vDotGradzi;
		wsub(t) = cv.budgetHeat.wSub( t);
		wsubAbove(t) = cv.budgetHeat.wSubAbove( t);
	end
	
	% TEMP: SAVE ADVECTION TERMS SEPARATELY
	%			cv.temp.meanU( t) = meanU;
	%			cv.temp.meanV( t) = meanV;
	%			cv.temp.dthetal_dx( t) = dthetal_dx;
	%			cv.temp.dthetal_dy( t) = dthetal_dy;
	%			cv.temp.dthetal_inv( t) = cv.avg.theta_l(  cv.avg.idxInversionBase(  t )+1, t) - cv.avg.theta_l(  cv.avg.idxInversionBase(  t )-1, t);
	%			cv.temp.dthetal_inv2( t) = cv.avg.theta_l(  cv.avg.idxInversionBase(  t )+1, t) - meanThetaL;
	%			cv.temp.meanThetaL( t) = meanThetaL;
	%			cv.temp.theta_l_inv( t) = cv.avg.theta_l(  cv.avg.idxInversionBase(  t )+1, t);
	%			cv.temp.dqt_dx( t) = dqt_dx;
	%			cv.temp.dqt_dy( t) = dqt_dy;
	%			cv.temp.dqt_inv( t) = cv.avg.qt(  cv.avg.idxInversionBase(  t )+1, t) - cv.avg.qt(  cv.avg.idxInversionBase(  t )-1, t);
	%			cv.temp.dqt_inv2( t) = cv.avg.qt(  cv.avg.idxInversionBase(  t )+1, t) - meanqt;
	%			cv.temp.meanqt( t) = meanqt;
	%			cv.temp.qt_inv( t) = cv.avg.qt(  cv.avg.idxInversionBase(  t )+1, t);
	
	% Advection of zi and ML-deviation term
	cv.avg.budgetMoisture.term1( t) = ( 1 ./ cv.avg.zi(  t ) ) .* vDotGradzi .* ( qt_zi_minus - meanqt ); % (1 / zi) .* vDotGradzi .* ( q_t_zi - q_t_BL )
	% Advection of q_t
	cv.avg.budgetMoisture.term2( t) = - vDotGradqt; % - vDotGradqt
	% Entrainment drying
% 	cv.avg.budgetMoisture.term3( t) = ( 1 ./ cv.avg.zi(  t ) ) .* w_e .* ( cv.avg.qt(  cv.avg.idxInversionBase(  t )+1, t) - cv.avg.qt(  cv.avg.idxInversionBase(  t )-1, t) ); % w_e / zi .* dqt_inv
	cv.avg.budgetMoisture.term3( t) = ( 1 ./ cv.avg.zi(  t ) ) .* w_e .* ( cv.avg.dqt( t ) ); % w_e / zi .* dqt_inv
	% Surface moistening
	% 			cv.budgetMoisture.term4( t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.lhf(  t ) ./ Lv ./ densityWater; %
	cv.avg.budgetMoisture.term4( t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.lhf(  t ) ./ Lv ./ 1.2; % HY TEMPORARY CHECK WITH DENSITY OF AIR RATHER THAN WATER
	% 			cv.budgetMoisture.term4(i-1, j-1, t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.qfx(  t ) ./ densityWater; % Use QFX rather than LHF (e.g. when using RUC)
	% Rain leaving PBL
	if t == 1
		cv.avg.budgetMoisture.term5( t) = - ( 1 ./ cv.avg.zi(  t ) ) .* (cv.avg.RainSFC(  t+1 ) - cv.avg.RainSFC(  t) ) ./ dt ./ 1000 .* densityWater ./ densityAir; % rain rate is negative
	else
		cv.avg.budgetMoisture.term5( t) = - ( 1 ./ cv.avg.zi(  t ) ) .* (cv.avg.RainSFC(  t ) - cv.avg.RainSFC(  t-1 )) ./ dt ./ 1000 .* densityWater ./ densityAir; % rain rate is negative
	end
	% Compute heating rates and moisture tendency using centered finite difference (forward and backward at boundaries)
	if t == 1
		cv.avg.budgetHeat.dTheta_L_dt_FD(  t ) = (1 ./ dt ) .* ( mean( cv.avg.theta_l(  1:cv.avg.idxInversionBase(  t ), t+1), 1 ) - mean( cv.avg.theta_l(  1:cv.avg.idxInversionBase(  t ), t ), 1 ) );
		cv.avg.budgetHeat.dTheta_dt_FD(  t ) = (1 ./ dt ) .* ( mean( cv.avg.theta(  1:cv.avg.idxInversionBase(  t ), t+1), 1 ) - mean( cv.avg.theta(  1:cv.avg.idxInversionBase(  t ), t ), 1 ) );
		cv.avg.budgetMoisture.dqt_dt_FD(  t ) = (1 ./ dt ) .* ( mean( cv.avg.qt(  1:cv.avg.idxInversionBase(  t ), t+1), 1 ) - mean( cv.avg.qt(  1:cv.avg.idxInversionBase(  t ), t ), 1 ) );
	elseif t == length(cv.time)
		cv.avg.budgetHeat.dTheta_L_dt_FD(  t ) = (1 ./ dt ) .* ( mean( cv.avg.theta_l(  1:cv.avg.idxInversionBase(  t ), t), 1 ) - mean( cv.avg.theta_l(  1:cv.avg.idxInversionBase(  t ), t-1 ), 1 ) );
		cv.avg.budgetHeat.dTheta_dt_FD(  t ) = (1 ./ dt ) .* ( mean( cv.avg.theta(  1:cv.avg.idxInversionBase(  t ), t), 1 ) - mean( cv.avg.theta(  1:cv.avg.idxInversionBase(  t ), t-1 ), 1 ) );
		cv.avg.budgetMoisture.dqt_dt_FD(  t ) = (1 ./ dt ) .* ( mean( cv.avg.qt(  1:cv.avg.idxInversionBase(  t ), t), 1 ) - mean( cv.avg.qt(  1:cv.avg.idxInversionBase(  t ), t-1 ), 1 ) );
	else
		cv.avg.budgetHeat.dTheta_L_dt_FD(  t ) = (1 ./ (2.*dt) ) .* ( mean( cv.avg.theta_l(  1:cv.avg.idxInversionBase(  t ), t+1), 1 ) - mean( cv.avg.theta_l(  1:cv.avg.idxInversionBase(  t ), t-1 ), 1 ) );
		cv.avg.budgetHeat.dTheta_dt_FD(  t ) = (1 ./ (2.*dt) ) .* ( mean( cv.avg.theta(  1:cv.avg.idxInversionBase(  t ), t+1), 1 ) - mean( cv.avg.theta(  1:cv.avg.idxInversionBase(  t ), t-1 ), 1 ) );
		cv.avg.budgetMoisture.dqt_dt_FD(  t ) = (1 ./ (2.*dt) ) .* ( mean( cv.avg.qt(  1:cv.avg.idxInversionBase(  t ), t+1), 1 ) - mean( cv.avg.qt(  1:cv.avg.idxInversionBase(  t ), t-1 ), 1 ) );
	end
	
	% Cloud base tendencies
	Tcb = cv.avg.T( cv.avg.idxCloudBase( t), t); % Temperature T at cloud base
	cv.avg.cloudBase.dzb_dqt( t) = (R .* Tcb) ./ (g .* meanqt) .* (1 - (Lv.*R)./(cpAir.*Rv.*Tcb)).^(-1); % Change in cloud base due to change in q_t (i.e. efficiency)
	cv.avg.cloudBase.dzb_dthetal( t) = 1 ./ ( g ./ (cpAir .* cv.avg.exner( cv.avg.idxCloudBase( t), t)) .* (1 - (cpAir.*Rv.*Tcb)./(R.*Lv)) ); % Change in cloud base due to change in theta_l (i.e. efficiency)
	% Compute cloud base tendency using centered finite difference (forward and backward at boundaries)
	if t == 1
		cv.avg.cloudBase.dzb_dt_FD( t) = (1 ./ dt ) .* ( cv.avg.zb( t+1) - cv.avg.zb( t) );
	elseif t == length(cv.time)
		cv.avg.cloudBase.dzb_dt_FD( t) = (1 ./ dt ) .* ( cv.avg.zb( t) - cv.avg.zb( t-1) );
	else
		cv.avg.cloudBase.dzb_dt_FD( t) = (1 ./ (2.*dt) ) .* ( cv.avg.zb( t+1) - cv.avg.zb( t-1) );
	end
end

%% Sum up to obtain heating rate and moisture tendency!
cv.avg.budgetHeat.dTheta_L_dt = cv.avg.budgetHeat.term1 + cv.avg.budgetHeat.term2 + cv.avg.budgetHeat.term3 + cv.avg.budgetHeat.term4 + cv.avg.budgetHeat.term5 + cv.avg.budgetHeat.term6;
cv.avg.budgetMoisture.dqt_dt = cv.avg.budgetMoisture.term1 + cv.avg.budgetMoisture.term2 + cv.avg.budgetMoisture.term3 + cv.avg.budgetMoisture.term4 + cv.avg.budgetMoisture.term5;

% Make assumption that residual between heating rate/moisture tendency is solely due to the entrainment term (i.e. all other terms are correct)
% --> back out entrainment term
cv.avg.budgetHeat.term3closure = cv.avg.budgetHeat.dTheta_L_dt_FD - cv.avg.budgetHeat.term1 - cv.avg.budgetHeat.term2 - cv.avg.budgetHeat.term4 - cv.avg.budgetHeat.term5 - cv.avg.budgetHeat.term6;
cv.avg.budgetMoisture.term3closure = cv.avg.budgetMoisture.dqt_dt_FD - cv.avg.budgetMoisture.term1 - cv.avg.budgetMoisture.term2 - cv.avg.budgetMoisture.term4 - cv.avg.budgetMoisture.term5;

%% Use FD HR to back out w_e from HB and MB residuals
% Use these to verify residual only contains entrainment term
for t = 1:length(cv.time)
% 			cv.avg.budgetMoisture.backedOutw_e( t) = ( cv.avg.budgetMoisture.dqt_dt_FD( t) - cv.avg.budgetMoisture.term1( t) - cv.avg.budgetMoisture.term2( t) - cv.avg.budgetMoisture.term4( t) - cv.avg.budgetMoisture.term5( t) ) .* cv.avg.zi(  t ) ./ ( cv.avg.qt(  cv.avg.idxInversionBase(  t )+1, t) - cv.avg.qt(  cv.avg.idxInversionBase(  t )-1, t) ); % dqt_dt .* z_i ./ dqt_inv
			cv.avg.budgetMoisture.backedOutw_e( t) = ( cv.avg.budgetMoisture.dqt_dt_FD( t) - cv.avg.budgetMoisture.term1( t) - cv.avg.budgetMoisture.term2( t) - cv.avg.budgetMoisture.term4( t) - cv.avg.budgetMoisture.term5( t) ) .* cv.avg.zi(  t ) ./ ( cv.avg.dqt(t) ); % dqt_dt .* z_i ./ dqt_inv
% 			cv.avg.budgetHeat.term3FromMoistureBudget( t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.budgetMoisture.backedOutw_e( t) .* ( cv.avg.theta_l(  cv.avg.idxInversionBase(  t )+1, t) - cv.avg.theta_l(  cv.avg.idxInversionBase(  t )-1, t) );
cv.avg.budgetHeat.term3FromMoistureBudget( t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.budgetMoisture.backedOutw_e( t) .* ( cv.avg.dthl(t) );
			
% 			cv.avg.budgetHeat.backedOutw_e( t) = ( cv.avg.budgetHeat.dTheta_L_dt_FD( t) - cv.avg.budgetHeat.term1( t) - cv.avg.budgetHeat.term2( t) - cv.avg.budgetHeat.term4( t) - cv.avg.budgetHeat.term5( t) - cv.avg.budgetHeat.term6( t) ) .* cv.avg.zi(  t ) ./ ( cv.avg.theta_l(  cv.avg.idxInversionBase(  t )+1, t) - cv.avg.theta_l(  cv.avg.idxInversionBase(  t )-1, t) );
			cv.avg.budgetHeat.backedOutw_e( t) = ( cv.avg.budgetHeat.dTheta_L_dt_FD( t) - cv.avg.budgetHeat.term1( t) - cv.avg.budgetHeat.term2( t) - cv.avg.budgetHeat.term4( t) - cv.avg.budgetHeat.term5( t) - cv.avg.budgetHeat.term6( t) ) .* cv.avg.zi(  t ) ./ ( cv.avg.dthl(t) );
% 			cv.avg.budgetHeat.term3FromHeatBudget( t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.budgetHeat.backedOutw_e( t) .* ( cv.avg.theta_l(  cv.avg.idxInversionBase(  t )+1, t) - cv.avg.theta_l(  cv.avg.idxInversionBase(  t )-1, t) );
			cv.avg.budgetHeat.term3FromHeatBudget( t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.budgetHeat.backedOutw_e( t) .* ( cv.avg.dthl(t) );
% 			cv.avg.budgetMoisture.term3FromHeatBudget( t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.budgetHeat.backedOutw_e( t) .* ( cv.avg.qt(  cv.avg.idxInversionBase(  t )+1, t) - cv.avg.qt(  cv.avg.idxInversionBase(  t )-1, t) ); % w_e / zi .* dqt_inv
			cv.avg.budgetMoisture.term3FromHeatBudget( t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.budgetHeat.backedOutw_e( t) .* ( cv.avg.dqt(t) ); % w_e / zi .* dqt_inv
			
			cv.avg.budgetHeat.meanw_e( t) = mean( [cv.avg.budgetHeat.backedOutw_e( t), cv.avg.budgetMoisture.backedOutw_e( t)] );
% 			cv.avg.budgetHeat.term3mean( t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.budgetHeat.meanw_e( t) .* ( cv.avg.theta_l(  cv.avg.idxInversionBase(  t )+1, t) - cv.avg.theta_l(  cv.avg.idxInversionBase(  t )-1, t) );
			cv.avg.budgetHeat.term3mean( t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.budgetHeat.meanw_e( t) .* ( cv.avg.dthl(t) );
% 			cv.avg.budgetMoisture.term3mean( t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.budgetHeat.meanw_e( t) .* ( cv.avg.qt(  cv.avg.idxInversionBase(  t )+1, t) - cv.avg.qt(  cv.avg.idxInversionBase(  t )-1, t) ); % w_e / zi .* dqt_inv
			cv.avg.budgetMoisture.term3mean( t) = ( 1 ./ cv.avg.zi(  t ) ) .* cv.avg.budgetHeat.meanw_e( t) .* ( cv.avg.dqt(t) ); % w_e / zi .* dqt_inv
end

% Sum up
cv.avg.budgetMoisture.dqt_dt_HB = cv.avg.budgetMoisture.term1 + cv.avg.budgetMoisture.term2 + cv.avg.budgetMoisture.term3FromHeatBudget + cv.avg.budgetMoisture.term4 + cv.avg.budgetMoisture.term5;
cv.avg.budgetMoisture.dqt_dt_mean = cv.avg.budgetMoisture.term1 + cv.avg.budgetMoisture.term2 + cv.avg.budgetMoisture.term3mean + cv.avg.budgetMoisture.term4 + cv.avg.budgetMoisture.term5;
cv.avg.budgetHeat.dTheta_L_dt_MB = cv.avg.budgetHeat.term1 + cv.avg.budgetHeat.term2 + cv.avg.budgetHeat.term3FromMoistureBudget + cv.avg.budgetHeat.term4 + cv.avg.budgetHeat.term5 + cv.avg.budgetHeat.term6;
cv.avg.budgetHeat.dTheta_L_dt_HB = cv.avg.budgetHeat.term1 + cv.avg.budgetHeat.term2 + cv.avg.budgetHeat.term3FromHeatBudget + cv.avg.budgetHeat.term4 + cv.avg.budgetHeat.term5 + cv.avg.budgetHeat.term6;
cv.avg.budgetHeat.dTheta_L_dt_mean = cv.avg.budgetHeat.term1 + cv.avg.budgetHeat.term2 + cv.avg.budgetHeat.term3mean + cv.avg.budgetHeat.term4 + cv.avg.budgetHeat.term5 + cv.avg.budgetHeat.term6;

%% Cloud base tendency
% dzb/dt = dzb/dthetal .* dthetal/dt + dzb/dqt .* dqt/dt;
cv.avg.cloudBase.dzb_dt_heat = cv.avg.cloudBase.dzb_dthetal .* cv.avg.budgetHeat.dTheta_L_dt_FD;
cv.avg.cloudBase.dzb_dt_moisture = cv.avg.cloudBase.dzb_dqt .* cv.avg.budgetMoisture.dqt_dt_FD;
cv.avg.cloudBase.dzb_dt = cv.avg.cloudBase.dzb_dt_heat + cv.avg.cloudBase.dzb_dt_moisture; % Use FD heating rate/moisture tendency

% Separate
cv.avg.cloudBase.dzb_dt_h1 = cv.avg.cloudBase.dzb_dthetal .* cv.avg.budgetHeat.term1;
cv.avg.cloudBase.dzb_dt_h2 = cv.avg.cloudBase.dzb_dthetal .* cv.avg.budgetHeat.term2;
cv.avg.cloudBase.dzb_dt_h3 = cv.avg.cloudBase.dzb_dthetal .* cv.avg.budgetHeat.term3mean;
cv.avg.cloudBase.dzb_dt_h4 = cv.avg.cloudBase.dzb_dthetal .* cv.avg.budgetHeat.term4;
cv.avg.cloudBase.dzb_dt_h5 = cv.avg.cloudBase.dzb_dthetal .* cv.avg.budgetHeat.term5;
cv.avg.cloudBase.dzb_dt_h6 = cv.avg.cloudBase.dzb_dthetal .* cv.avg.budgetHeat.term6;

cv.avg.cloudBase.dzb_dt_m1 = cv.avg.cloudBase.dzb_dqt .* cv.avg.budgetMoisture.term1;
cv.avg.cloudBase.dzb_dt_m2 = cv.avg.cloudBase.dzb_dqt .* cv.avg.budgetMoisture.term2;
cv.avg.cloudBase.dzb_dt_m3 = cv.avg.cloudBase.dzb_dqt .* cv.avg.budgetMoisture.term3mean;
cv.avg.cloudBase.dzb_dt_m4 = cv.avg.cloudBase.dzb_dqt .* cv.avg.budgetMoisture.term4;
cv.avg.cloudBase.dzb_dt_m5 = cv.avg.cloudBase.dzb_dqt .* cv.avg.budgetMoisture.term5;

for i = 2:(size(cv.avg.u, 1)-1)
	for j = 2:(size(cv.avg.u, 2)-1)
		cv.avg.cloudBase.wrfThickness( 1) = cv.avg.zi( 1) - cv.avg.zb( 1);
		cv.avg.cloudBase.predictedThickness( 1) = cv.avg.cloudBase.wrfThickness( 1);
		cv.avg.cloudBase.predictedzb( 1) = cv.avg.zb( 2);
		cv.avg.cloudBase.predictedzi( 1) = cv.avg.zi( 2);
		for t = 2:length(cv.time)
			cv.avg.cloudBase.predictedzi( t) = cv.avg.cloudBase.predictedzi( t-1) + (cv.avg.budgetHeat.meanw_e( t-1) - cv.avg.budgetHeat.wSub( t-1)) .* dt; % dzi/dt + adv = w_e + w_sub
			cv.avg.cloudBase.predictedzb( t) = cv.avg.cloudBase.predictedzb( t-1) + cv.avg.cloudBase.dzb_dt( t-1) .* dt;
			if cv.avg.cloudBase.predictedzb( t) < 0, cv.avg.cloudBase.predictedzb( t) = 0; end
			cv.avg.cloudBase.predictedThickness( t) = cv.avg.zi( t) - cv.avg.cloudBase.predictedzb( t);
			if cv.avg.cloudBase.predictedThickness( t) < 0, cv.avg.cloudBase.predictedThickness( t) = 0; end
			cv.avg.cloudBase.wrfThickness( t) = cv.avg.zi( t) - cv.avg.zb( t);
			if cv.avg.cloudBase.wrfThickness( t) < 0, cv.avg.cloudBase.wrfThickness( t) = 0; end
		end
		
	end
end

%%%%%%%% END FUNCTION

% End function
fprintf('ML CV done.\n')

%% Write MLM soundings if requested
if makeSounding
	
	% sound_in
	% height (m)	pressure (Pa)	Theta_l (K)	qt (g/kg)
	fid = fopen([soundingDir '/sound_in_' datestr(startTime, 'yyyymmdd')], 'w');
	fprintf(fid, '0.0 %f %f %f\n', cv.spatial.psfc(i, j, 1), cv.spatial.tl2(i, j, 1), 1000.*cv.spatial.q2(i, j, 1));
	for z = 1:size(cv.heightAGL, 3)
		fprintf(fid, '%f %f %f %f\n', cv.heightAGL(i, j, z, 1), cv.spatial.pressure(i, j, z, 1), cv.spatial.theta_l(i, j, z, 1), 1000.*cv.spatial.qt(i, j, z, 1));
	end
	fclose(fid);
	
	% advection_in
	% NOTE TO SELF:	Experiment with doing dT/dx as (meanT1 - meanT2) and mean(T1vec - T2vec). Hopefully, they will be the same value, or close. If not, (meanT1 - meanT2) should be more
	%				in line with mixed-layer derivation.
	% time	Theta_l_advection (K/s)	qt_advection (g/kg/s) zi_advection (m/s)
	
	% Process time into doy format [PST]
	doyTime = doy( cv.time );
	
	fid = fopen([soundingDir '/advection_in_' datestr(startTime, 'yyyymmdd')], 'w');
	for t = 1:length( doyTime )
		fprintf(fid, '%f %.20f %.20f %.20f\n', doyTime(t), thetaLAdvection(t), qtAdvection(t), ziAdvection(t));
	end
	fclose(fid);
	
	fid = fopen([soundingDir '/subsidence_in_' datestr(startTime, 'yyyymmdd')], 'w');
	for t = 1:length( doyTime )
		fprintf(fid, '%f %.20f %.20f %.20f\n', doyTime(t), wsub(t), wsubAbove(t));
	end
	fclose(fid);
	
	% NAMELIST
	fid = fopen([soundingDir '/NAMELIST_' datestr(startTime, 'yyyymmdd')], 'w');
	fprintf(fid, 'strtim %f\n', doyTime(1));
	fprintf(fid, 'timmax 61200\n'); % Max time until we run out of advection data
	fprintf(fid, 'cntlat %f\n', wrf.XLAT(x(i), y(j), 1));
	fprintf(fid, 'cntlon %f\n', wrf.XLONG(x(i), y(j), 1));
	fprintf(fid, 'theta0 289\n');
	fprintf(fid, 'iratyp 2\n');
	fprintf(fid, 'subtyp 2\n');
	fprintf(fid, 'advtyp 2\n');
	fprintf(fid, 'wmatyp 1\n');
	fprintf(fid, 'isftyp 1\n');
	fprintf(fid, 'bowenr %f\n', mean( squeeze([cv.spatial.shf(i, j, :)]) ./ squeeze([cv.spatial.lhf(i, j, :)]) ) );
	fprintf(fid, 'srfeff 0.88\n');
	fprintf(fid, 'Csknav 20000\n');
	fprintf(fid, 'tsrfce %f\n', cv.spatial.tsk(i, j, 1));
	fclose(fid);
end

end