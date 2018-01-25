%% Run after go.m to check convective velocity scales
% Formulate after Fang et al. (2013) "Turbulence in Continental Stratocumulus, Part I: External Forcings and Turbulence Structures"

% W_t* ^3 = W_s*^3 + W_r*^3, where
% W_s* = [ (g*z_t / theta_v_bar) * (w'theta_v')_bar ] ^ 1/3, where (w'theta_v')_bar is surface buoyancy flux
% W_r* = [ (g*z_t * (- dF_r)) / (theta_v_bar * rho * cp) ] ^ 1/3, where - dF_r is the difference in total net radiative flux between cloud top and cloud base

n = length(ts.avg.wstar);
cp = 1006; % [J / kg K]

for idx = 1:n
	wstar_sfc(idx) = ts.avg.wstar_sfc{idx};
	wstar_deardoff(idx) = ts.avg.wstar{idx};
	cloudBase(idx) = ts.avg.zbmn{idx};
	cloudTop(idx) = ts.avg.zcmn{idx};
end

% Convention: F_net positive towards surface
for idx = 1:n
	topIdx = findMatch( ps.zt, cloudTop(idx) );
	botIdx = findMatch( ps.zt, cloudBase(idx) );
	if ~isempty(topIdx) && ~isempty(botIdx)
		fNetTop(idx) = ps.avg.sflxd{idx}(topIdx) - ps.avg.sflxu{idx}(topIdx) + ps.avg.lflxd{idx}(topIdx) - ps.avg.lflxu{idx}(topIdx);
		fNetBot(idx) = ps.avg.sflxd{idx}(botIdx) - ps.avg.sflxu{idx}(botIdx) + ps.avg.lflxd{idx}(botIdx) - ps.avg.lflxu{idx}(botIdx);
	else
		fNetTop(idx) = 0; fNetBot(idx) = 0;
	end
	wstar_radiative(idx) = ( (9.8 .* ts.avg.(conf.zScale){idx} .* ( - (fNetTop(idx) - fNetBot(idx)) ) ) ./ (ts.avg.tsrf{idx} .* rho .* cp) ) .^ (1/3);
end

%% Plots
% % Check F_net,top and F_net,bot
% figure; plot(fNetTop), hold on, plot(fNetBot, 'b--')
% legend('F_{net,top}', 'F_{net,bot}')

% Total net radiative flux between cloud top and cloud base
figure; plot(fNetTop - fNetBot), legend('F_{net,top} - F_{net,bot}')
ylabel('Heat Flux [W m^{-2}]')
xlabel('Time (hours)')

% Velocity scales
figure; plot(real( wstar_sfc + wstar_radiative ) , 'k')
hold on % ; plot(real(wstar_deardoff), 'r')
plot(real(wstar_sfc))
plot(real(wstar_radiative), 'm')
plot(cloudBase./1000, 'g--')
plot(cloudTop./1000, 'g')
ylabel('Velocity [m s^{-1}]')
xlabel('Time (hours)')
legend('w*_{tot}', 'w*_{sfc}', 'w*_{rad}', 'z_b' ,'z_t')

% % Rate of change of velocity scales
% figure; plot(diff(real(wstar_sfc)))
% hold on; plot(diff(real(wstar_deardoff)), 'r')
% plot(cloudBase./1000, 'g')
% plot(cloudTop./1000, 'g--')