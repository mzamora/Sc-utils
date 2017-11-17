function [p z tp td RH w winddir windspeed theta]=getSounding(date_i)
% [p z tp td RH w winddir windspeed theta theta_e theta_v]=getSounding(date_i)
% Get sounding variables below 3km, skipping mandatory levels. 

    
    %% Load yearly netcdf file (obtained from https://ruc.noaa.gov/raobs/ )
    raob_filename=['raw/raob_soundings_',num2str(year(date_i)),'.cdf'];
    raob_dates=datetime(1970,1,1,0,0,0)+ncread(raob_filename,'synTime')/3600/24;
    i_star=find(raob_dates==date_i);
    
    %% Read all necessary variables (cut below 3km, since we are interested in low level properties)
    htSigT=ncread(raob_filename,'htSigT',[1 i_star],[250 1]); htSigT(htSigT>1e30)=[]; f=htSigT<3000; htSigT=htSigT(f);
    tpSigT=ncread(raob_filename,'tpSigT',[1 i_star],[250 1]); tpSigT(tpSigT>1e30)=[]; tpSigT=tpSigT(f); tpSigT(tpSigT==99999)=nan;
    tdSigT=ncread(raob_filename,'tdSigT',[1 i_star],[250 1]); tdSigT(tdSigT>1e30)=[]; tdSigT=tdSigT(f); tdSigT(tdSigT==99999)=nan;
    wdSigT=ncread(raob_filename,'wdSigT',[1 i_star],[250 1]); wdSigT(wdSigT>1e30)=[]; wdSigT=wdSigT(f);
    wsSigT=ncread(raob_filename,'wsSigT',[1 i_star],[250 1]); wsSigT(wsSigT>1e30)=[]; wsSigT=wsSigT(f);
    prSigT=ncread(raob_filename,'prSigT',[1 i_star],[250 1]); prSigT(prSigT>1e30)=[]; prSigT=prSigT(f);
    
    htMan=ncread(raob_filename,'htMan',[1 i_star],[22 1]); htMan(htMan>1e30)=[]; f=htMan<3000; htMan=htMan(f);
    tpMan=ncread(raob_filename,'tpMan',[1 i_star],[22 1]); tpMan(tpMan>1e30)=[]; tpMan=tpMan(f); tpMan(tpMan==99999)=nan;
    tdMan=ncread(raob_filename,'tdMan',[1 i_star],[22 1]); tdMan(tdMan>1e30)=[]; tdMan=tdMan(f); tdMan(tdMan==99999)=nan;
    wdMan=ncread(raob_filename,'wdMan',[1 i_star],[22 1]); wdMan(wdMan>1e30)=[]; wdMan=wdMan(f);
    wsMan=ncread(raob_filename,'wsMan',[1 i_star],[22 1]); wsMan(wsMan>1e30)=[]; wsMan=wsMan(f);
    prMan=ncread(raob_filename,'prMan',[1 i_star],[22 1]); prMan(prMan>1e30)=[]; prMan=prMan(f);
    
    htSigW=ncread(raob_filename,'htSigW',[1 i_star],[150 1]); htSigW(htSigW>1e30)=[]; f=htSigW<3000; htSigW=htSigW(f);
    tpSigW=ncread(raob_filename,'tpSigW',[1 i_star],[150 1]); tpSigW(tpSigW>1e30)=[]; tpSigW=tpSigW(f); tpSigW(tpSigW==99999)=nan;
    tdSigW=ncread(raob_filename,'tdSigW',[1 i_star],[150 1]); tdSigW(tdSigW>1e30)=[]; tdSigW=tdSigW(f); tdSigW(tdSigW==99999)=nan;
    wdSigW=ncread(raob_filename,'wdSigW',[1 i_star],[150 1]); wdSigW(wdSigW>1e30)=[]; wdSigW=wdSigW(f);
    wsSigW=ncread(raob_filename,'wsSigW',[1 i_star],[150 1]); wsSigW(wsSigW>1e30)=[]; wsSigW=wsSigW(f);
    prSigW=ncread(raob_filename,'prSigW',[1 i_star],[150 1]); prSigW(prSigW>1e30)=[]; prSigW=prSigW(f);
    
    %% Combine SigT, SigW and the surface point in Man to create the final set
    ht=[htMan(1); htSigT]; 
    tp=[tpMan(1); tpSigT];
    td=[tdMan(1); tdSigT];
    wd=[wdMan(1); wdSigT];
    ws=[wsMan(1); wsSigT];
    pr=[prMan(1); prSigT];
    
    %% Sort all based on height
    [z,ies]=sort(ht);
    tp=tp(ies); td=td(ies); wd=wd(ies); ws=ws(ies); pr=pr(ies);
    
    %% Get extra variables
    exner=(pr/1e5).^0.286;
    theta=tp./exner;
    TP=tp-273.15; TD=TP-td; %temps in celsius
    es=6.1094*exp(17.625*TP./(TP+243.04));
    e=6.1094*exp(17.625*TD./(TD+243.04));
    RH=e./es;
    w=0.622*e./(pr-e);
    wsat=0.622*es./(pr-es);
    
end