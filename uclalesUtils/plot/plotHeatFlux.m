function plotHeatFlux(ps, segment, z, zi_bar, selections, lw, fs, style)

hold on
set(gca, 'fontsize', fs.axes)

count = 1;
for idx = selections
	plot(ps.avg.tot_tw{idx} ./ ps.avg.tot_tw{idx}(1), z./zi_bar{idx}, style{count}, 'linewidth', lw)
	timestamps{count} = sprintf('t = %.1f hr', segment.etime(idx)./3600);
	count = count + 1;
end

% Plot resolved and SGS for last time step
plot(ps.avg.sfs_tw{ selections(end) } ./ ps.avg.sfs_tw{ selections(end) }(1), z./zi_bar{ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '-.', 'linewidth', lw);
plot( ( ps.avg.tot_tw{ selections(end) } - ps.avg.sfs_tw{ selections(end) } ) ./ ps.avg.sfs_tw{ selections(end) }(1), z./zi_bar{ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '--', 'linewidth', lw);

axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('<w''\Theta''>/<w_{0}''\Theta_{0}''> [-]', 'fontsize', fs.axes)
box on
legend([timestamps, 'SGS (last)', 'Res (last)'], 'location', 'northeastoutside', 'fontsize', fs.legend)

end