function plotBuoyProd(ps, ts, conf, selections, lw, fs, style)

hold on
set(gca, 'fontsize', fs.axes)

count = 1;
for idx = selections

	BuoyProd = ps.avg.boy_prd{idx}; % Buoyancy production

	% Plot
	plot(BuoyProd, ps.zm./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw)
	timestamps{count} = sprintf('t = %.1f hr', ps.segment.etime(idx)./3600);
	count = count + 1;
end

SGS = (ps.avg.sfs_boy{idx}); % Normalize by something?

% Plot resolved and SGS for last time step
plot(SGS, ps.zm./ts.avg.(conf.zScale){ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '-.', 'linewidth', lw);
plot( (BuoyProd - SGS) , ps.zm./ts.avg.(conf.zScale){ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '--', 'linewidth', lw);

axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes)
xlabel('Buoyancy production of resolved TKE [m^{2}/s^{3}]', 'fontsize', fs.axes)
box on
legend([timestamps, 'SGS_{end}', 'Res_{end}'], 'location', 'northeastoutside', 'fontsize', fs.legend)

end