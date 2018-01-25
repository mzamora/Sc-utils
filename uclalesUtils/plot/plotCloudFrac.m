function plotCloudFrac(ts, lw, fs)

set(gca, 'fontsize', fs.axes)

plot(ts.time./3600, ts.cfrac, 'k', 'linewidth', lw)

axis tight
ylabel('Cloud fraction [-]', 'fontsize', fs.axes)
xlabel('Time [hr]', 'fontsize', fs.axes)
box on
ylim([0 1])
end