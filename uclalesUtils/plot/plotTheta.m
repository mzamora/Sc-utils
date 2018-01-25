function plotTheta(theta, segment, z, zi_bar, selections, lw, fs, style)

hold on
set(gca, 'fontsize', fs.axes)

count = 1;
for idx = selections
	plot(theta{idx}, z./zi_bar{idx}, style{count}, 'linewidth', lw)
	timestamps{count} = sprintf('t = %.1f hr', segment.etime(idx)./3600);
	count = count + 1;
end

axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('\Theta [K]', 'fontsize', fs.axes)
box on

legend(timestamps, 'location', 'best', 'fontsize', fs.legend)

end