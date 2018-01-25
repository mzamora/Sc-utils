function plotSurfBuoyancy(ts, lw, fs)

set(gca, 'fontsize', fs.axes)

plot(ts.time./3600, ts.sfcbflx, 'k', 'linewidth', lw)

axis tight
ylabel('Surface buoyancy flux [m/s^{2}]', 'fontsize', fs.axes)
xlabel('Time [hr]', 'fontsize', fs.axes)
box on
end