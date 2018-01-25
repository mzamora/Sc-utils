function plotTKE(time, tke, lw, fs)

set(gca, 'fontsize', fs.axes)

plot(time./3600, tke, 'k', 'linewidth', lw)

axis tight
ylabel('Integrated TKE [m^3/s^2]', 'fontsize', fs.axes)
xlabel('Time [hr]', 'fontsize', fs.axes)
box on
end