function plotFluxW(ps, ts, conf, selections, lw, fs, style)

hold on
set(gca, 'fontsize', fs.axes)

count = 1;
for idx = selections

	FluxW = (ps.avg.tot_ww{idx})./ ts.avg.(conf.uScale){idx}.^2; % <w'w'> ./ w*^2

	% Plot
	plot(FluxW, ps.zm./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw)
	timestamps{count} = sprintf('t = %.1f hr', ps.segment.etime(idx)./3600);
	count = count + 1;
end

SGS = (ps.avg.sfs_ww{idx})./ ts.avg.(conf.uScale){idx}.^2; % <w'w'> ./ w*^2

% Plot resolved and SGS for last time step
plot(SGS, ps.zm./ts.avg.(conf.zScale){ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '-.', 'linewidth', lw);
plot( (FluxW - SGS) , ps.zm./ts.avg.(conf.zScale){ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '--', 'linewidth', lw);

axis tight
% ylim([0.8 1.4])
ylim([0 2])
xlim([0 0.6])
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes)
xlabel('<w''w''>/w^{*}^{2} [-]', 'fontsize', fs.axes)
box on
legend([timestamps, 'SGS_{end}', 'Res_{end}'], 'location', 'northeastoutside', 'fontsize', fs.legend)

end