function plotLWP(ts, lw, fs)

set(gca, 'fontsize', fs.axes)

hold on
plot(ts.time./3600, ts.lwp_bar, 'k', 'linewidth', lw)

axis tight
ylabel('Mean liquid water path [g/m^{2}]', 'fontsize', fs.axes)
xlabel('Time [hr]', 'fontsize', fs.axes)
box on
end