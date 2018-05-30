function plotOblRi(ps, ts, conf, lw, fs, style, h)
defineConstants

for idx = 1:length(ps.avg.t)
    vecUstar(idx) = ts.avg.ustar{idx};
    vecTsfc(idx) = ps.avg.theta_v{idx}(1);
    idxZi = findMatch(ps.zt, ts.avg.zi1_bar{idx});
    vecTpbl(idx) = ps.avg.theta_v{idx}(idxZi - 10);
    vecDthv(idx) = ps.avg.theta_v{idx}(idxZi+1) - ps.avg.theta_v{idx}(idxZi-1);
    vecDz(idx) = ps.zt(idxZi+1) - ps.zt(idxZi-1);
    vecDU(idx) = ps.avg.u{idx}(idxZi+1) - ps.avg.u{idx}(idxZi-1);
    vecDV(idx) = ps.avg.v{idx}(idxZi+1) - ps.avg.v{idx}(idxZi-1);
    vecWTV(idx) = ts.avg.sfcbflx{idx} ./ 9.81 .* vecTsfc(idx);
end

Obl = - (vecUstar .^ 3 .* vecTsfc) ./ (0.4 .* 9.81 .* vecWTV);
Rib = ( (9.81 ./ vecTpbl) .* (vecDthv) .* (vecDz) ) ./ (vecDU.^2 + vecDV.^2);
% Rig = ;


% Generate plots
figure(h(1))
hold on, set(gca, 'fontsize', fs.axes)
plot(ps.segment.stime./3600, Obl)
figure(h(2))
hold on, set(gca, 'fontsize', fs.axes)
plot(ps.segment.stime./3600, Rib)


% Plots
figure(h(1))
axis tight
% ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('<u''^{2}>/w*^{2} [-]', 'fontsize', fs.axes)
ylabel('Obl [m]'), xlabel('Time [h]')
box on, grid on
% legend(timestamps, 'location', 'northeast', 'fontsize', fs.legend)
% set(gca, 'ytick', [0:0.2:2]), ylim([0.8 1.4])
% xlim([0 0.8])

figure(h(2))
axis tight
% ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes), xlabel('<v''^{2}>/w*^{2} [-]', 'fontsize', fs.axes)
ylabel('Ri [-]'), xlabel('Time [h]')
box on, grid on
% legend(timestamps, 'location', 'northeast', 'fontsize', fs.legend)
% set(gca, 'ytick', [0:0.2:2]), ylim([0.8 1.4])
% xlim([0 0.8])

end