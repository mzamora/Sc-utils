function plotVarProfiles(ps, ts, conf, selections, lw, fs, style, h)

count = 1;

for idx = selections
	figure(h(1))
	hold on, set(gca, 'fontsize', fs.axes)
	plot(ps.avg.u_2{idx} ./ ts.avg.(conf.uScale){idx}.^2, ps.zt./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw)
	figure(h(2))
	hold on, set(gca, 'fontsize', fs.axes)
	plot(ps.avg.v_2{idx} ./ ts.avg.(conf.uScale){idx}.^2, ps.zt./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw)
	figure(h(3))
	hold on, set(gca, 'fontsize', fs.axes)
	plot(ps.avg.w_2{idx} ./ ts.avg.(conf.uScale){idx}.^2, ps.zt./ts.avg.(conf.zScale){idx}, style{count}, 'linewidth', lw)
	timestamps{count} = sprintf('t = %.1f hr', ps.segment.etime(idx)./3600);
	count = count + 1;
end

figure(h(1))
axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('<u''^{2}>/w*^{2} [-]', 'fontsize', fs.axes)
box on, grid on
legend(timestamps, 'location', 'northeast', 'fontsize', fs.legend)
set(gca, 'ytick', [0:0.2:2]), ylim([0.8 1.4])
xlim([0 0.8])

figure(h(2))
axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('<v''^{2}>/w*^{2} [-]', 'fontsize', fs.axes)
box on, grid on
legend(timestamps, 'location', 'northeast', 'fontsize', fs.legend)
set(gca, 'ytick', [0:0.2:2]), ylim([0.8 1.4])
xlim([0 0.8])

figure(h(3))
axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('<w''^{2}>/w*^{2} [-]', 'fontsize', fs.axes)
box on, grid on
legend(timestamps, 'location', 'northeast', 'fontsize', fs.legend)
set(gca, 'ytick', [0:0.2:2]), ylim([0.8 1.4])
xlim([0 0.4])

end