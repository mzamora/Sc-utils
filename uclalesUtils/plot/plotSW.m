function plotSW(ts, lw, fs)

set(gca, 'fontsize', fs.axes)

hold on
plot(ts.time./3600, ts.sflxut, 'k', 'linewidth', lw)
plot(ts.time./3600, ts.sflxdt, '--k', 'linewidth', lw)

plot(ts.time./3600, ts.sflxus, 'r', 'linewidth', lw)
plot(ts.time./3600, ts.sflxds, '--r', 'linewidth', lw)

axis tight
ylabel('Shortwave radiative flux [W/m^{2}]', 'fontsize', fs.axes)
xlabel('Time [hr]', 'fontsize', fs.axes)
legend('TOA SW up', 'TOA SW down', 'Sfc SW up', 'Sfc SW down')
box on
end