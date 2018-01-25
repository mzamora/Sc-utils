% Radiative flux vs. height (total, shortwave, longwave [subtract short from total])

function plotRadFlux(ps, ts, conf, selections, lw, fs, style)

hold on
set(gca, 'fontsize', fs.axes)

count = 1;
for idx = selections

	TotalRad = ps.avg.rflx{idx}; % Total radiative flux

	% Plot
	plot(TotalRad, ps.zm./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw)
	timestamps{count} = sprintf('t = %.1f hr', ps.segment.etime(idx)./3600);
	count = count + 1;
end

SWFlux = (ps.avg.sflx{idx}); % Shortwave flux

% Plot SW and LW for last time step
plot(SWFlux, ps.zm./ts.avg.(conf.zScale){ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '-.', 'linewidth', lw);
plot( (TotalRad - SWFlux) , ps.zm./ts.avg.(conf.zScale){ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '--', 'linewidth', lw);

axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes)
xlabel('Radiative flux [W/m^{2}]', 'fontsize', fs.axes)
box on
legend([timestamps, 'SW_{end}', 'LW_{end}'], 'location', 'northeastoutside', 'fontsize', fs.legend)

end