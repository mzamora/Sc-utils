function plotScalarVarProfiles(ps, ts, conf, selections, lw, fs, style, h)
defineConstants;
count = 1;

for idx = selections
	figure(h(1))
	hold on, set(gca, 'fontsize', fs.axes)
% 	plot(ps.avg.t_2{idx} ./ (max(ps.avg.tot_tw{idx}) ./ (cpAir .* densityAir .* ts.avg.(conf.uScale){idx}) ).^2, ps.zt./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw)
	plot(ps.avg.t_2{idx}, ps.zt./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw)
	
	figure(h(2))
	hold on, set(gca, 'fontsize', fs.axes)
% 	plot(ps.avg.q_2{idx} ./ (max(ps.avg.tot_qw{idx}) ./ Lv ./ densityAir ./ ts.avg.(conf.uScale){idx} ).^2, ps.zt./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw)
	plot(ps.avg.q_2{idx}, ps.zt./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw)
	
	timestamps{count} = sprintf('t = %.1f hr', ps.segment.etime(idx)./3600);
	
	figure(h(3))
	yyaxis left
	hold on, set(gca, 'fontsize', fs.axes)
	plot(idx, trapz(ps.zt, ps.avg.t_2{idx}), '*')
	
	yyaxis right
	plot(idx, trapz(ps.zt, ps.avg.q_2{idx}), '+')
	count = count + 1;
end

figure(h(1))
axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('<\theta''^{2}> [K^2]', 'fontsize', fs.axes)
box on, grid on
legend(timestamps, 'location', 'northeast', 'fontsize', fs.legend)
set(gca, 'ytick', [0:0.2:2]), ylim([0.8 1.4])
xlim([0 4])

figure(h(2))
axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('<q_t''^{2}> [kg/kg]', 'fontsize', fs.axes)
box on, grid on
legend(timestamps, 'location', 'northeast', 'fontsize', fs.legend)
set(gca, 'ytick', [0:0.2:2]), ylim([0.8 1.4])
xlim([0 3e-6])

figure(h(3))
yyaxis left
ylabel('I_\theta [K^2 m]')
xlabel('Time [hr]')
xlim([0.9*selections(1) 1.1*selections(end)])
% ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('<q_t''^{2}> [kg/kg]', 'fontsize', fs.axes)
yyaxis right
ylabel('I_{qt} [kg m/ kg]')
xlim([0.9*selections(1) 1.1*selections(end)])
title('* = I_\theta, + = I_{qt}')
box on, grid on
% legend('I_\theta', 'I_{qt}', 'location', 'northeast', 'fontsize', fs.legend)
% set(gca, 'ytick', [0:0.2:2]), ylim([0.8 1.4])


end