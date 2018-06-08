function plotFluxQ(ps, segment, z, zi_bar, selections, lw, fs, style)

hold on
set(gca, 'fontsize', fs.axes)

count = 1;
for idx = selections
	plot(ps.avg.tot_qw{idx} ./ ps.avg.tot_qw{idx}(1), z./zi_bar{idx}, style{count}, 'linewidth', lw)
	timestamps{count} = sprintf('t = %.1f hr', segment.etime(idx)./3600);
	count = count + 1;
end

% Plot resolved and SGS for last time step
plot(ps.avg.sfs_qw{ selections(end) } ./ ps.avg.sfs_qw{ selections(end) }(1), z./zi_bar{ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '-.', 'linewidth', lw);
plot( ( ps.avg.tot_qw{ selections(end) } - ps.avg.sfs_qw{ selections(end) } ) ./ ps.avg.sfs_qw{ selections(end) }(1), z./zi_bar{ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '--', 'linewidth', lw);

axis tight
ylim([0.8 1.4])
% ylim([0 2])
xlim([0 4])
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('<w''q''>/<w_{0}''q_0''> [-]', 'fontsize', fs.axes)
box on
legend([timestamps, 'SGS (last)', 'Res (last)'], 'location', 'northeastoutside', 'fontsize', fs.legend)

end