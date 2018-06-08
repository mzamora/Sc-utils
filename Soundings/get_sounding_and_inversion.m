function [p z tp td RH w winddir windspeed theta theta_e theta_v zi ni zit nit]=get_sounding_and_inversion(date_i,scutilsdir)
% [p z tp td RH w winddir windspeed theta theta_e theta_v zi ni zit nit]=get_sounding_and_inversion(date_i)
% Returns complete sounding variables and inversion base and top
% Monica Zamora, UCSD SRAF (Dec 2017) http://solar.ucsd.edu

    %% Get all sounding variables
    pathName=[scutilsdir,'/Soundings/raw/72293_',datestr(date_i,'yyyy_mm_'),'0100_',num2str(eomday(year(date_i),month(date_i))),'12.csv']; %set source of data
    ListofVar={'PRES';'HGHT';'TEMP';'DWPT';'RELH';'MIXR';'WDIR';'WSPD';'THTA';'THTE';'THTV'};
    [p z tp td RH w winddir windspeed theta theta_e theta_v]=Get_sounding_Var(date_i,date_i+0.5,pathName,ListofVar);
    %% Get inversion
    [~,~,zit,zi,nit,ni]=TMP_Inversion_Strength_Cal(tp,z/1000,z(1)); %Xiaohui's code for inversion height
    [~,~,hght_top2,~,eta_top2,~]=TMP_Inversion_Strength_Cal(-w,z/1000,z(1)); %Xiaohui's code for inversion height
    if hght_top2<zit
        zit=hght_top2;
        nit=eta_top2;
    end
            
end