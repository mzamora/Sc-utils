%% How does the thermodynamic fit for qt, tl, zi perform?

% Read LES file data
LESdata=csvread('/mnt/Sc_group/LES_analysis/timeSeries/NKX_20170712_ts_BL_thetaL_qt_zi_zb.csv',1,0);

f=~LESdata(:,2)==0; %filter out times with no data
time=LESdata(f,1)/3600; %time in hh
tl_BL=LESdata(f,2); 
qt_BL=LESdata(f,3);
zi=LESdata(f,4);
zb=LESdata(f,5);

% Create fit data
a=182.1; b=0.0041; c=0.672; % fit coefficients
qt_fit=-a-b*zi+c*tl_BL;
tl_fit=(b*zi+qt_BL+a)/c;
zi_fit=(c*tl_BL-qt_BL-a)/b;

% Compute cloud base from zi fitted
zb_derived=get_wellmixed_zb(zi_fit,qt_BL,tl_BL);

% Plot fits
subplot(511); plot(time,zi,time,zi_fit); legend('LES','Fit z_i(q_t,θ_l)'); xlabel('Time [hh]'); ylabel('z_i [m]')
subplot(512); plot(time,qt_BL,time,qt_fit); legend('LES','Fit q_t(z_i,θ_l)'); xlabel('Time [hh]'); ylabel('q_t^{BL} [g/kg]')
subplot(513); plot(time,tl_BL,time,tl_fit); legend('LES','Fit θ_l(z_i,q_t)'); xlabel('Time [hh]'); ylabel('θ_l^{BL} [K]')
subplot(514); plot(time,zb,time,zb_derived); legend('LES','From z_i fit'); xlabel('Time [hh]'); ylabel('z_b [m]')
subplot(515); plot(time,zi-zb,time,zi_fit-zb_derived'); legend('LES','From z_i fit'); xlabel('Time [hh]'); ylabel('h [m]')
