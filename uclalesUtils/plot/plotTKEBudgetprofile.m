function plotTKEBudgetprofile(ps, ts, conf, selections, lw, fs, style)

hold on
set(gca, 'fontsize', fs.axes)

if (length(selections) > 1); error('Only choose 1 selection.'); end

normalization = ts.avg.(conf.uScale){ selections }.^3 ./ ts.avg.(conf.zScale){ selections};

% Shear
plot( ps.avg.shr_prd{ selections } + ps.avg.sfs_shr{ selections }, ps.zt./ts.avg.(conf.zScale){selections}, 'g', 'linewidth', lw);
plot( ps.avg.sfs_shr{ selections }, ps.zt./ts.avg.(conf.zScale){selections}, 'g--', 'linewidth', lw); % SFS
% Buoy
plot( ps.avg.boy_prd{ selections } + ps.avg.sfs_boy{ selections }, ps.zt./ts.avg.(conf.zScale){selections}, 'color', [.98 .60 .01], 'linewidth', lw);
plot( ps.avg.sfs_boy{ selections }, ps.zt./ts.avg.(conf.zScale){selections}, 'color', [.98 .6 .01], 'linestyle', '--', 'linewidth', lw); % SFS
% Transport
plot( ps.avg.trans{ selections }, ps.zt./ts.avg.(conf.zScale){selections}, 'b', 'linewidth', lw);
% Dissipation
plot( - ps.avg.diss{ selections }, ps.zt./ts.avg.(conf.zScale){selections}, 'r', 'linewidth', lw);

axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('TKE budget / (w*^3/z_i^3) [-]', 'fontsize', fs.axes)
box on
grid on
legend({'Shear_{r+s}', 'Shear_{sgs}', 'Buoy_{r+s}', 'Buoy_{sgs}', 'Trans', 'Diss'}, 'location', 'northeastoutside', 'fontsize', fs.legend)
xlim([-0.002 0.004])

end