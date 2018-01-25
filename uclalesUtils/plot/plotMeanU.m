% Shortwave radiative flux up and down

function plotMeanU(ps, ts, conf, selections, lw, fs, style)

hold on
set(gca, 'fontsize', fs.axes)

count = 1;
for idx = selections
	meanU = ps.avg.u{idx}; % LW up
	
	% Plot
	plot(meanU./ts.avg.ustar{idx}, ps.zm./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw) % <u>/u*
	timestamps{count} = sprintf('t = %.1f hr', ps.segment.etime(idx)./3600);
	count = count + 1;
end

axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes)
xlabel('<u>/u* [-]', 'fontsize', fs.axes)
box on
legend(timestamps, 'location', 'northeastoutside', 'fontsize', fs.legend)

end