function plotLW(ts, lw, fs)

set(gca, 'fontsize', fs.axes)

hold on
plot(ts.time./3600, ts.lflxut, 'k', 'linewidth', lw)
plot(ts.time./3600, ts.lflxdt, '--k', 'linewidth', lw)

plot(ts.time./3600, ts.lflxus, 'r', 'linewidth', lw)
plot(ts.time./3600, ts.lflxds, '--r', 'linewidth', lw)

axis tight
ylabel('Longwave radiative flux [W/m^{2}]', 'fontsize', fs.axes)
xlabel('Time [hr]', 'fontsize', fs.axes)
legend('TOA LW up', 'TOA LW down', 'Sfc LW up', 'Sfc LW down')
box on
end