% To read variables in NAM (from panthers)
run /home/mzamoraz/nctoolbox-1.1.3/setup_nctoolbox.m

%% Describe site properties
cases={'Gilman'}; 
lats=32.877394; %latitudes
lons=-117.233874; %longitudes
NAMfolder='/mnt/lab_48tb2d/database/NAM_ARCHIVE/'; %from panthers

for icase=1:length(cases)
    % get coordinates close to each location
    [iis(icase),jjs(icase)]=find_NAM_indices(lats(icase),lons(icase));
end

%% Process wanted years: will save 1 file per each year
for yy=2015:2018

%initialize variables
iday=1; failed_times=[];
ghi=nan(372,24,length(cases));
RH=nan(372,24,length(cases));
windu=nan(372,24,length(cases));
windv=nan(372,24,length(cases));
temp_2m=nan(372,24,length(cases));

for day=datetime(yy-1,12,30+(0:371)) %4 extra days for time difference or averaging
    for ih=0:23 %24h data
        ih6=floor(ih/6)*6; % 00, 06, 12, 18 utc folders only
        ihfile=ih-ih6;
        if day>datetime(2017,2,23) %filename change in the archive UwU
            NAMfile=[NAMfolder,num2str(ih6,'%02i'),'_UTC/',datestr(day,'YYYY/mm/'),datestr(day,'yyyymmdd'),'_',num2str(ihfile,'%02i'),'.tm00.grib2'];
            else
            NAMfile=[NAMfolder,num2str(ih6,'%02i'),'_UTC/',datestr(day,'YYYY/mm/'),datestr(day,'yyyymmdd'),'_',num2str(ihfile,'%02i'),'.grb2.tm00'];
        end
        try %load data
            nc = ncgeodataset(NAMfile);
            %% read data % to read variables: nc.variables
            for icase=1:length(cases)
                temp_2m(iday,ih+1,icase)=nc.data('Temperature_height_above_ground',[1,1,iis(icase),jjs(icase)],[1,1,iis(icase),jjs(icase)]);
                RH(iday,ih+1,icase)=nc.data('Relative_humidity_height_above_ground',[1,1,iis(icase),jjs(icase)],[1,1,iis(icase),jjs(icase)]);
                ghi(iday,ih+1,icase)=nc.data('Downward_Short-Wave_Radiation_Flux_surface',[1,iis(icase),jjs(icase)],[1,iis(icase),jjs(icase)]);
                windu(iday,ih+1,icase)=nc.data('u-component_of_wind_height_above_ground',[1,1,iis(icase),jjs(icase)],[1,1,iis(icase),jjs(icase)]);
                windv(iday,ih+1,icase)=nc.data('v-component_of_wind_height_above_ground',[1,1,iis(icase),jjs(icase)],[1,1,iis(icase),jjs(icase)]);
            end % cases
            clear nc            
        catch %errors
            fprintf(['Error day ',num2str(iday),' h ',num2str(ih),' \n'])
            failed_times=[failed_times,day+ih/24];
%            exist(NAMfile)
        end
    end
        iday=iday+1

%% Save 1h data
save(['Gilman_NAM',num2str(yy),'.mat'],'temp_srf','temp_2m','ghi','RH','failed_times')
end %day
end %year

%% %%%%%%%%%%%%% START FROM HERE IF ONLY POST PROCESSING %%%%%%%%%%%%%%%%%%
%% Write 5 min data, appending to the GHI data
load('../data/Temperature_NAM.mat');
temp_2mfilled=temp_2m; % to fill nans
for icase=1: length(cases)
    % missing days
    temp_sum=sum(temp_2m(:,:,icase),2); 
    missing=find(isnan(temp_sum));
    missingblocks=find(diff(missing)~=1); i0=1;
    for iblock=1:length(missingblocks) % treat each block separately
        block=temp_2m(missing(i0)-2:missing(missingblocks(iblock))+2,:,icase); %extended block to interpolate
        temp_t=reshape(block',1,size(block,1)*size(block,2)); 
        fnan0=find(isnan(temp_t),1)-1; fnan1=find(isnan(temp_t),1,'last')+1;

        for ih=1:24 %we interpolate each hour column
            block(:,ih)=fillmissing(block(:,ih),'spline');
        end
        % smooth the block (only the nan part)
        temp_t=reshape(block',1,size(block,1)*size(block,2)); temp_t(fnan0:fnan1)=smooth(temp_t(fnan0:fnan1),'sgolay');
        block=reshape(temp_t,size(block,2),size(block,1))';     
        % put the block back
        temp_2mfilled(missing(i0)-2:missing(missingblocks(iblock))+2,:,icase)=block;
        i0=missingblocks(iblock)+1; %start of next block
    end
    
    % arrange 1h temp and time
    temp=reshape(temp_2mfilled(:,:,icase)',1,372*24)-273.15; % to Celsius
    temp_time=datetime(2017,12,28,(0:372*24-1),0,0)-timediff(icase)/24;
    
    %% 5 min data
    temp_f=temp(temp_time>=datetime(2018,1,1) & temp_time<datetime(2019,1,1));
    time_5min=datetime(2018,1,1,0,0,0):minutes(5):datetime(2018,12,31,23,55,0);
    temp_5min=resample(temp_f,12,1);
        
    load(['../data/Results_ensuredmean',strrep(cases{icase},' ','_')])
    %plot(temp_time,temp,Results.Timestamp,Results.Power/200) %check plot
    Results.Temperature_C=temp_5min;
    save(['../data/Results_ensuredmean',strrep(cases{icase},' ','_')],'Results','Data','SAData')
    
    %% write updated table
    T=table(Results.Timestamp',Results.Power',Results.Temperature_C','VariableNames',{'Time','AdjustedPower_kW','Temperature_C'});
    writetable(T,['../analysis/Results_ensuredmean',strrep(cases{icase},' ','_'),'.csv'],'Delimiter',',')
end

