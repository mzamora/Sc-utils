function plotSurfHeatFlux(ts, lw, fs)

set(gca, 'fontsize', fs.axes)

hold on
plot(ts.time./3600, ts.shf_bar, 'k', 'linewidth', lw)
plot(ts.time./3600, ts.lhf_bar, '--k', 'linewidth', lw)

axis tight
ylabel('Surface heat flux [W/m^{2}]', 'fontsize', fs.axes)
xlabel('Time [hr]', 'fontsize', fs.axes)
legend('Sensible', 'Latent', 'location', 'southeast')
box on
end