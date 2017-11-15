function [p z tp td RH w winddir windspeed theta theta_e theta_v]=getSounding(date_i)
    
    %% Load yearly netcdf file (obtained from https://ruc.noaa.gov/raobs/ )
    raob_filename=['raw/raob_soundings_',num2str(year(date_i)),'.cdf'];
    raob_dates=datetime(1970,1,1,0,0,0)+ncread(raob_filename,'synTime')/3600/24;
    i_star=find(raob_dates==date_i);
    
    %% Read all necessary variables
    htSigT=ncread(filename,'htSigT',[1 i_star],[250 1]); htSigT(htSigT>1e30)=[]; f=htSigT<3000; htSigT=htSigT(f);
    tpSigT=ncread(filename,'tpSigT',[1 i_star],[250 1]); tpSigT(tpSigT>1e30)=[]; tpSigT=tpSigT(f); tpSigT(tpSigT==99999)=nan;
    tdSigT=ncread(filename,'tdSigT',[1 i_star],[250 1]); tdSigT(tdSigT>1e30)=[]; tdSigT=tdSigT(f); tdSigT(tdSigT==99999)=nan;
    wdSigT=ncread(filename,'wdSigT',[1 i_star],[250 1]); wdSigT(wdSigT>1e30)=[]; wdSigT=wdSigT(f);
    wsSigT=ncread(filename,'wsSigT',[1 i_star],[250 1]); wsSigT(wsSigT>1e30)=[]; wsSigT=wsSigT(f);
    prSigT=ncread(filename,'prSigT',[1 i_star],[250 1]); prSigT(prSigT>1e30)=[]; prSigT=prSigT(f);
    
    htMan=ncread(filename,'htMan',[1 i_star],[22 1]); htMan(htMan>1e30)=[]; f=htMan<3000; htMan=htMan(f);
    tpMan=ncread(filename,'tpMan',[1 i_star],[22 1]); tpMan(tpMan>1e30)=[]; tpMan=tpMan(f); tpMan(tpMan==99999)=nan;
    tdMan=ncread(filename,'tdMan',[1 i_star],[22 1]); tdMan(tdMan>1e30)=[]; tdMan=tdMan(f); tdMan(tdMan==99999)=nan;
    wdMan=ncread(filename,'wdMan',[1 i_star],[22 1]); wdMan(wdMan>1e30)=[]; wdMan=wdMan(f);
    wsMan=ncread(filename,'wsMan',[1 i_star],[22 1]); wsMan(wsMan>1e30)=[]; wsMan=wsMan(f);
    prMan=ncread(filename,'prMan',[1 i_star],[22 1]); prMan(prMan>1e30)=[]; prMan=prMan(f);
    
    htSigW=ncread(filename,'htSigW',[1 i_star],[150 1]); htSigW(htSigW>1e30)=[]; f=htSigW<3000; htSigW=htSigW(f);
    tpSigW=ncread(filename,'tpSigW',[1 i_star],[150 1]); tpSigW(tpSigW>1e30)=[]; tpSigW=tpSigW(f); tpSigW(tpSigW==99999)=nan;
    tdSigW=ncread(filename,'tdSigW',[1 i_star],[150 1]); tdSigW(tdSigW>1e30)=[]; tdSigW=tdSigW(f); tdSigW(tdSigW==99999)=nan;
    wdSigW=ncread(filename,'wdSigW',[1 i_star],[150 1]); wdSigW(wdSigW>1e30)=[]; wdSigW=wdSigW(f);
    wsSigW=ncread(filename,'wsSigW',[1 i_star],[150 1]); wsSigW(wsSigW>1e30)=[]; wsSigW=wsSigW(f);
    prSigW=ncread(filename,'prSigW',[1 i_star],[150 1]); prSigW(prSigW>1e30)=[]; prSigW=prSigW(f);
    
    
    
    
    
    
end