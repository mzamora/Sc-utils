function plotTKEprofile(ps, ts, conf, selections, lw, fs, style)

hold on
set(gca, 'fontsize', fs.axes)

count = 1;
for idx = selections
	TKE = 0.5 * (ps.avg.u_2{idx} + ps.avg.v_2{idx} + ps.avg.w_2{idx});
	plot(TKE./ts.avg.(conf.uScale){idx}.^2, ps.zt./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw)
	timestamps{count} = sprintf('t = %.1f hr', ps.segment.etime(idx)./3600);
	count = count + 1;
end

% Plot resolved and SGS for last time step
plot( ( TKE - ps.avg.sfs_tke{ selections(end) } ) ./ ts.avg.(conf.uScale){ selections(end) }.^2, ps.zt./ts.avg.(conf.zScale){ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '--', 'linewidth', lw); % resolved
plot(ps.avg.sfs_tke{ selections(end) } ./ ts.avg.(conf.uScale){ selections(end) }.^2, ps.zt./ts.avg.(conf.zScale){ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '-.', 'linewidth', lw); % SFS


axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('TKE/w*^2 [-]', 'fontsize', fs.axes)
box on
legend([timestamps, 'Res (last)', 'SGS (last)'], 'location', 'northeastoutside', 'fontsize', fs.legend)
ylim([0 2])
xlim([0 0.8])

end