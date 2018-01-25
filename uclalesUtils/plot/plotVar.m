function plotVar(ps, lw, fs)

hold on
set(gca, 'fontsize', fs.axes)

for idx = 1:size(ps.u_2, 2)
	plot(ps.time./3600, mean(ps.u_2,1), 'k', 'linewidth', lw)
	plot(ps.time./3600, mean(ps.v_2,1), 'r', 'linewidth', lw)
	plot(ps.time./3600, mean(ps.w_2,1), 'g', 'linewidth', lw)
end

axis tight
ylabel('Velocity variance [m^2/s^2]', 'fontsize', fs.axes), xlabel('Time [hr]', 'fontsize', fs.axes)
legend([{'u'}, {'v'}, {'w'}], 'location', 'northwest', 'fontsize', fs.legend)
box on

end