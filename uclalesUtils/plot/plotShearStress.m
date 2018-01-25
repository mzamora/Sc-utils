function plotShearStress(ps, ts, segment, z, zi_bar, selections, velScale, lw, fs, style, vecAvg)

hold on
set(gca, 'fontsize', fs.axes)

if nargin < 11
	vecAvg = 0;
end

count = 1;
for idx = selections
	% (<u’w’>^2 + <v’w’>^2)^¼ / u*_sfc
	% Compute vector average for resolved scale
% 	vectorAvgShearStress = ( ( (ps.avg.tot_uw{idx}).^2 + (ps.avg.tot_vw{idx}).^2 ) .^ (1/4) ) ./ ts.avg.(velScale){idx}; % sqrt(vecAvg) ./ w*
	if vecAvg == 0,
		ShearStress = (ps.avg.tot_uw{idx})./ ts.avg.(velScale){idx}.^2; % <u'w'> ./ w*^2
	else
		ShearStress = ( ( (ps.avg.tot_uw{idx}).^2 + (ps.avg.tot_vw{idx}).^2 ) .^ (1/2) ) ./ ts.avg.(velScale){idx}.^2; % vecAvg ./ w*^2
	end

	% Plot
	plot(ShearStress, z./zi_bar{idx}, style{count}, 'linewidth', lw)
	timestamps{count} = sprintf('t = %.1f hr', segment.etime(idx)./3600);
	count = count + 1;
end

% SGS = ( ( (ps.avg.sfs_uw{ selections(end) }).^2 + (ps.avg.sfs_vw{ selections(end) }).^2 ) .^ (1/4) ) ./ ts.avg.(velScale){ selections(end) }; % sqrt(vecAvg) ./ w*
if vecAvg == 0,
	SGS = (ps.avg.sfs_uw{idx})./ ts.avg.(velScale){idx}.^2; % <u'w'> ./ w*^2
else
	SGS = ( ( (ps.avg.sfs_uw{ selections(end) }).^2 + (ps.avg.sfs_vw{ selections(end) }).^2 ) .^ (1/2) ) ./ ts.avg.(velScale){ selections(end) }.^2; % vecAvg ./ w*^2
end

% Plot resolved and SGS for last time step
plot(SGS, z./zi_bar{ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '-.', 'linewidth', lw);
plot( (ShearStress - SGS) , z./zi_bar{ selections(end) }, 'color', [0.7 0.7 0.7], 'linestyle', '--', 'linewidth', lw);

axis tight
ylabel('z/<z_{i}> [-]', 'fontsize', fs.axes)

% if strcmpi(velScale, 'ustar'), xlabel('<u''w''>^{1/2}/u^{*} [-]', 'fontsize', fs.axes), elseif strcmpi(velScale, 'wstar'), xlabel('(<u''w''>^{2} + <v''w''>^2)^{1/4}/w^{*} [-]', 'fontsize', fs.axes), end
if vecAvg == 0,
	if strcmpi(velScale, 'ustar'), xlabel('<u''w''>/u^{*}^{2} [-]', 'fontsize', fs.axes), elseif strcmpi(velScale, 'wstar'), xlabel('<u''w''>/w^{*}^{2} [-]', 'fontsize', fs.axes), end
else
	if strcmpi(velScale, 'ustar'), xlabel('(<u''w''>^{2} + <v''w''>^2)^{1/2}/u^{*}^{2} [-]', 'fontsize', fs.axes), elseif strcmpi(velScale, 'wstar'), xlabel('(<u''w''>^{2} + <v''w''>^2)^{1/2}/w^{*}^{2} [-]', 'fontsize', fs.axes), end
end
box on
legend([timestamps, 'SGS_{end}', 'Res_{end}'], 'location', 'northeastoutside', 'fontsize', fs.legend)

end