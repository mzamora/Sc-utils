function plotWater(ps, ts, conf, selections, lw, fs, style)

hold on
set(gca, 'fontsize', fs.axes)

count = 1;
for idx = selections
	plot(ps.avg.q{idx}, ps.zt./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw)
	timestamps{count} = sprintf('t = %.1f hr', ps.segment.etime(idx)./3600);
	count = count + 1;
end

count = 1;
for idx = selections % Another for loop so legend is correct
	plot(ps.avg.l{idx}, ps.zt./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw, 'linestyle', '--')
	count = count + 1;
end

axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('Mixing ratio [g / kg]', 'fontsize', fs.axes)
box on
legend([timestamps, 'Liquid'], 'location', 'best', 'fontsize', fs.legend)

end