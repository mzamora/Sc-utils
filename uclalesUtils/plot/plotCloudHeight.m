function plotCloudHeight(ts, lw, fs)

set(gca, 'fontsize', fs.axes)

hold on
plot(ts.time./3600, ts.zcmn, 'k', 'linewidth', lw)
plot(ts.time./3600, ts.zbmn, '--k', 'linewidth', lw)

axis tight
ylabel('Mean cloud height [m]', 'fontsize', fs.axes)
xlabel('Time [hr]', 'fontsize', fs.axes)
legend('Top', 'Base', 'location', 'southeast')
ylim([0 max(ts.zcmn)+200])
box on
end