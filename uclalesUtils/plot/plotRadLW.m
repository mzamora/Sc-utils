% Longwave radiative flux up and down

function plotRadLW(ps, ts, conf, selections, lw, fs, style)

hold on
set(gca, 'fontsize', fs.axes)

count = 1;
for idx = selections
	LW_UP = ps.avg.lflxu{idx}; % LW up

	% Plot
	plot(LW_UP, ps.zm./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw)
	timestamps{count} = sprintf('t = %.1f hr', ps.segment.etime(idx)./3600);
	count = count + 1;
end

count = 1;
for idx = selections % Another for loop so legend is correct
	LW_DOWN = ps.avg.lflxd{idx}; % LW down
	plot(LW_DOWN, ps.zm./ts.avg.(conf.zScale){ selections(end) }, style{count}, 'linestyle', '--', 'linewidth', lw);
	count = count + 1;
end

axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes)
xlabel('Radiative flux [W/m^{2}]', 'fontsize', fs.axes)
box on
legend([timestamps, 'Down'], 'location', 'northeastoutside', 'fontsize', fs.legend)

end