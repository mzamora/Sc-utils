ind=1;
Cp=1004; Rd=287.05; Lv=2.5e6;

for yy=2014:2017
    for mm=5:9
        switch mm
            case {5,7,8}
                dmax=31;
            case {6,9}
                dmax=30;
        end
    
        for dd=1:dmax
            for hh=[0 12]
                pathName=strcat('raw/72293_',num2str(yy),'_',sprintf('%.2d',mm),'_0100_',sprintf('%.2d',dmax),'12.csv');
                VarList={'PRES';'HGHT';'TEMP';'DWPT';'RELH';'MIXR';'WDIR';'WSPD';'THTA';'THTE';'THTV'}; %list of parameters available
                date0=datetime(yy,mm,dd,hh,0,0);
                [p z tp td RH w winddir windspeed theta theta_e theta_v]=Get_sounding_Var(date0,date0+0.5,pathName,VarList); %Xiaohui's code
                p=double(p); %in hPa
                z=double(z);
                %%
                if isempty(z)
                    continue
                end
                
                %%
                date_i(ind,1)=date0;

                %% Filter sounding data
                zup=3000; % up to 3km
                f=isnan(w)|(z>zup)|(w>50); %Wipe out NaN values and data above desired height or weird values
                p(f)=[]; td(f)=[]; tp(f)=[]; w(f)=[]; theta(f)=[]; theta_e(f)=[]; theta_v(f)=[]; windspeed(f)=[]; winddir(f)=[]; z(f)=[]; RH(f)=[];
                
%                 figure(1); clf
%                 annotation('textbox', [0 0.9 1 0.1],'String',strcat('i=',num2str(ind),', ',datestr(date_i(ind))),'EdgeColor','none','HorizontalAlignment','center','FontSize',12)
%                 subplot(131); plot(theta,z,'-*'); hold on
%                 subplot(132); plot(w,z,'-*'); hold on
%                 subplot(133); plot(RH,z,'-*',[95 95],[z(1) z(end)],'r--'); hold on
                
                %% Find LCL
                try 
                    [~,~,LCL_bolton]=Cal_LCL_Bolton(tp,p,w,[],z,td);
%                     subplot(131); plot([min(theta) max(theta)],[LCL_bolton(1) LCL_bolton(1)],'r:')
%                     subplot(132); plot([min(w) max(w)],[LCL_bolton(1) LCL_bolton(1)],'r:')
%                     subplot(133); plot([min(RH) max(RH)],[LCL_bolton(1) LCL_bolton(1)],'r:')
                    LCL_srf(ind,1)=LCL_bolton(1);
                catch
                    LCL_srf(ind,1)=nan;
                end
                
                %% Cloud base and number of cloud segments from RH profile (last cloud)
                try
                    cloudpoints=find(RH>95);
                    n_clouds(ind,1)=1;
                    zb_index=cloudpoints(1);
                    z_cloudbase(ind,1)=z(cloudpoints(1));
                    for ic=1:length(cloudpoints)-1
                        if cloudpoints(ic+1)~=(cloudpoints(ic)+1)
                            n_clouds(ind,1)=n_clouds(ind)+1;
                            z_cloudbase(ind,1)=z(cloudpoints(ic+1));
                            zb_index=cloudpoints(ic+1);
                        end
                    end    
                    dLCL_decoupling(ind,1)=z_cloudbase(ind)-LCL_bolton(1);
                    LCLCB=LCL_bolton(find(LCL_bolton<z_cloudbase(ind),1,'last'));
                    if isempty(LCLCB)
                        LCL_CB(ind,1)=nan;
                    else
                        LCL_CB(ind,1)=LCLCB;
                    end
                    dtv1(ind,1)=theta_v(zb_index)-theta_v(1);
                    dtv2(ind,1)=theta_v(zb_index)-theta_v(2);
                catch
                    n_clouds(ind,1)=nan; z_cloudbase(ind,1)=nan; dLCL_decoupling(ind,1)=nan; LCL_CB(ind,1)=nan; zb_index=nan;
                    dtv1(ind,1)=nan; dtv2(ind,1)=nan;
                end
                
                %% Detect if there's fog present
                try
                    if(RH(1)>95)
                        fog(ind,1)=1;
                    else
                        fog(ind,1)=1;
                    end
                catch
                    fog(ind,1)=nan;
                end
                
                
                %% Find inversion
                try
                    [DT_max,numofInv,hght_top,hght_base,~,~]=TMP_Inversion_Strength_Cal(tp,z/1000,z(1)); 
                    [~,~,inv_w_top,inv_w_base,~,~]=TMP_Inversion_Strength_Cal(-w,z/1000,z(1)); %check inversion levels for mixing ratio profile
                    assert(~isnan(hght_base))
                    
                    z_inv_base(ind,1)=hght_base*1000; %Inversion base height (IBH)
                    z_inv_top(ind,1)=min(hght_top,inv_w_top)*1000; %Inversion top height
                    zi_index1=find(z==round(z_inv_base(ind)),1); %find index for IBH
                    if(zi_index1==1)
                        fprintf('Inversion at sfc')
                        sfc_inversion(ind,1)=1;
                        top_inversion(ind,1)=0;
                    elseif (z_inv_top(ind)==z(end)) %skip if IBH is at surface or at last point
                        fprintf('Inversion at top')
                        top_inversion(ind,1)=1;
                        sfc_inversion(ind,1)=1;
                    else
                        sfc_inversion(ind,1)=0;
                        top_inversion(ind,1)=0;
                    end
                    zi_index2=find(z==round(z_inv_top(ind),0));
%                     subplot(131); plot([min(theta) max(theta)],[z_inv_base(ind) z_inv_base(ind)],'r',[min(theta) max(theta)],[z_inv_top(ind) z_inv_top(ind)],'r--'); 
%                     subplot(132); plot([min(w) max(w)],[z_inv_base(ind) z_inv_base(ind)],'r',[min(w) max(w)],[z_inv_top(ind) z_inv_top(ind)],'r--'); 
%                     subplot(133); plot([min(RH) max(RH)],[z_inv_base(ind) z_inv_base(ind)],'r',[min(RH) max(RH)],[z_inv_top(ind) z_inv_top(ind)],'r');
                catch
                    z_inv_base(ind,1)=nan; z_inv_top(ind,1)=nan; sfc_inversion(ind,1)=nan; top_inversion(ind,1)=nan; zi_index1=nan; zi_index2=nan; 
                end
                
                %% Idealized well-mixed profile
                try
                    qT_BL(ind,1)=trapz(z(1:zi_index1),w(1:zi_index1))/(z_inv_base(ind)-z(1)); %Weighted average for the mixed layer
                    qT_3km(ind,1)=trapz(z(zi_index2:end),w(zi_index2:end))/(z(end)-z_inv_top(ind)); %Weighted average for the mixed layer
                    qT_jump(ind,1)=qT_3km(ind)-qT_BL(ind);
                    
                    if (~isnan(zb_index))&&(zb_index>1)
                        thetaL_BL(ind,1)=trapz(z(1:zb_index+1),theta(1:zb_index+1))/(z(zb_index+1)-z(1)); %Weighted average below cloud base
                    else
                        thetaL_BL(ind,1)=trapz(z(1:zi_index1),theta(1:zi_index1))/(z_inv_base(ind)-z(1)); %Weighted average for the mixed layer
                    end
                    poly=polyfit(z(zi_index2:end),theta(zi_index2:end),1); %Polynomial fit for above
                    thetaL_3km(ind,1)=polyval(poly,3000); %Alternative: original data theta_l(zi_index+1:end);
                    thetaL_jump(ind,1)=polyval(poly,z_inv_base(ind))-thetaL_BL(ind);

%                     subplot(131); plot([thetaL_BL(ind) thetaL_BL(ind)],[z(1) z_inv_base(ind)],'r',[thetaL_BL(ind)+thetaL_jump(ind) thetaL_3km(ind)],[z_inv_base(ind) 3000],'r'); 
%                     subplot(132); plot([qT_BL(ind) qT_BL(ind)],[z(1) z_inv_base(ind)],'r',[qT_3km(ind) qT_3km(ind)],[z_inv_base(ind) 3000],'r')
                catch
                    qT_BL(ind,1)=nan; qT_3km(ind,1)=nan; qT_jump(ind,1)=nan;
                    thetaL_BL(ind,1)=nan; thetaL_3km(ind,1)=nan; thetaL_jump(ind,1)=nan;
                end
                
%                 subplot(131); ylabel('Height [m]'); xlabel('Potential temperature [K]'); xlim([min(theta) max(theta)+1]); ylim([0 3000])
%                 subplot(132); xlabel('Mixing ratio [g kg^{-1}]'); ylabel('Height'); xlim([min(w) max(w)+.01]); ylim([0 3000])
%                 subplot(133); xlabel('RH [%]'); ylabel('Height'); xlim([0 100]); ylim([0 3000])
                
                %% Delta values for decoupling as in Jones
                if (~isnan(zi_index1))&&(zi_index1>5)
                    try
                    n2=round(zi_index1/3,0);
                    q_bot(ind,1)=trapz(z(1:n2),w(1:n2))/(z(n2)-z(1));
                    q_up(ind,1)=trapz(z(zi_index1-n2:zi_index1),w(zi_index1-n2:zi_index1))/(z_inv_base(ind)-z(zi_index1-n2));
                    theta_bot(ind,1)=trapz(z(1:n2),theta(1:n2))/(z(n2)-z(1));
                    theta_up(ind,1)=trapz(z(zi_index1-n2:zi_index1),theta(zi_index1-n2:zi_index1))/(z_inv_base(ind)-z(zi_index1-n2));
                    dq_decoupling(ind,1)=q_bot(ind)-q_up(ind);
                    dtheta_decoupling(ind,1)=theta_bot(ind)-theta_up(ind);
                    catch
                    q_bot(ind,1)=nan; q_up(ind,1)=nan; theta_bot(ind,1)=nan; theta_up(ind,1)=nan; dq_decoupling(ind,1)=nan; dtheta_decoupling(ind,1)=nan;
                    end
                else
                    q_bot(ind,1)=nan; q_up(ind,1)=nan; theta_bot(ind,1)=nan; theta_up(ind,1)=nan; dq_decoupling(ind,1)=nan; dtheta_decoupling(ind,1)=nan;
                end
                
                %% BL wind
                try
                    BLwnd_spd_avg(ind,1)=trapz(z(1:zi_index1),windspeed(1:zi_index1))/(z(zi_index1)-z(1));
                    u=windspeed.*sin(winddir*pi/180); BLu=trapz(z(1:zi_index1),u(1:zi_index1))/(z(zi_index1)-z(1));
                    v=windspeed.*cos(winddir*pi/180); BLv=trapz(z(1:zi_index1),v(1:zi_index1))/(z(zi_index1)-z(1));
                    BLwnd_dir_avg(ind,1)=mod(atan(BLu/BLv)*180/pi,360);
                catch
                    BLwnd_spd_avg(ind,1)=nan; BLwnd_dir_avg(ind,1)=nan;
                end
                %%
%                 text(60,2700,strcat('n_{cld}=',num2str(n_clouds(ind))))
%                 text(60,2500,strcat('Δtv=',num2str(dtv1(ind))))
%                 text(60,2300,strcat('Δzb=',num2str(dLCL_decoupling(ind))))
%                 text(60,2100,strcat('Δq=',num2str(dq_decoupling(ind))))
%                 text(60,1900,strcat('Δθ=',num2str(dtheta_decoupling(ind))))
%                 print(strcat('allsoundings/',datestr(date_i(ind),'yyyymmdd-hh')),'-depsc')
                ind=ind+1;

            end
        end
    end
end

%% Create table
wellmixed_dtv=(dtv1<.25);
decoupled_dtv=(dtv1>1);
wellmixed_dLCL=dLCL_decoupling<150;
decoupled_dLCL=dLCL_decoupling>150;

ALL_SOUNDINGS=table(date_i,n_clouds,LCL_srf,z_cloudbase,dLCL_decoupling,z_inv_base, ...
    z_inv_top,sfc_inversion,top_inversion,qT_BL,qT_jump,qT_3km,thetaL_BL, ...
    thetaL_jump,thetaL_3km,BLwnd_dir_avg,BLwnd_spd_avg,dq_decoupling,dtheta_decoupling,wellmixed_dtv, ...
    decoupled_dtv,wellmixed_dLCL,decoupled_dLCL);

ALL_SOUNDINGS.Properties.Description='Sounding parameters for NKX, May-Sept 2014-2017';
ALL_SOUNDINGS.Properties.VariableUnits={'UTC time' '' 'm' 'm' 'm' 'm' 'm' '' '' 'g/kg' 'g/kg' 'g/kg' 'K' 'K' 'K' '°' 'm/s' 'g/kg' 'K' 'K' 'K' 'm' 'm' };
ALL_SOUNDINGS.Properties.VariableDescriptions{1}='Sounding Time';
ALL_SOUNDINGS.Properties.VariableDescriptions{2}='Number of cloud layer from the RH profile (RH>95%). Cases with 0 clouds are not wanted';
ALL_SOUNDINGS.Properties.VariableDescriptions{3}='LCL computed using the Bolton formula';
ALL_SOUNDINGS.Properties.VariableDescriptions{4}='Cloud base of the last layer of cloud based on the RH profile';
ALL_SOUNDINGS.Properties.VariableDescriptions{5}='Difference between LCL and the cloud base. Measure of decoupling following Jones';
ALL_SOUNDINGS.Properties.VariableDescriptions{6}='Inversion base height';
ALL_SOUNDINGS.Properties.VariableDescriptions{7}='Inversion top height';
ALL_SOUNDINGS.Properties.VariableDescriptions{8}='1 if the inversion found is at the surface (not wanted)';
ALL_SOUNDINGS.Properties.VariableDescriptions{9}='1 if the inversion found is at the top of the profile (not wanted)';
ALL_SOUNDINGS.Properties.VariableDescriptions{10}='Total water mixing ratio in the boundary layer. Well mixed assumption';
ALL_SOUNDINGS.Properties.VariableDescriptions{11}='Total water mixing ratio jump at the top of the boundary layer. Well mixed assumption';
ALL_SOUNDINGS.Properties.VariableDescriptions{12}='Total water mixing ratio at 3km. Well mixed assumption';
ALL_SOUNDINGS.Properties.VariableDescriptions{13}='Liquid water potential temperature in the boundary layer. Well mixed assumption';
ALL_SOUNDINGS.Properties.VariableDescriptions{14}='Liquid water potential temperature jump at the top of the boundary layer. Well mixed assumption';
ALL_SOUNDINGS.Properties.VariableDescriptions{15}='Liquid water potential temperature a 3km. Well mixed assumption';
ALL_SOUNDINGS.Properties.VariableDescriptions{16}='Direction of the boundary layer wind velocity average';
ALL_SOUNDINGS.Properties.VariableDescriptions{17}='Speed of the boundary layer wind velocity average';
ALL_SOUNDINGS.Properties.VariableDescriptions{18}='Δq between bottom and top of the BL. Measure for decoupling following Ghate';
ALL_SOUNDINGS.Properties.VariableDescriptions{19}='Δθv between bottom and top of the BL. Measure for decoupling following Ghate';
ALL_SOUNDINGS.Properties.VariableDescriptions{20}='1: it is probably well-mixed since Δθv<0.25';
ALL_SOUNDINGS.Properties.VariableDescriptions{21}='1: it is probably decoupled since Δθv>1';
ALL_SOUNDINGS.Properties.VariableDescriptions{22}='1: it is probably well-mixed since ΔLCL<150';
ALL_SOUNDINGS.Properties.VariableDescriptions{23}='1: it is probably decoupled since ΔLCL>150';

save('NKX_soundings_table','ALL_SOUNDINGS')
writetable(ALL_SOUNDINGS,'NKX_soundings_table.csv')