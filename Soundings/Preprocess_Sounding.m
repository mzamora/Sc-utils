%This function is used to preprocess Sounding data.

%Remeber to do human checks to improve the quality of data because there
%are always missing data.

% addpath(genpath('/home/xiaohui/WRF/m_programe_Linux'));
% addpath(genpath('/home/xiaohui/REST2/m_programe'));

%The dirName is the folder containing sounding data.
dirName=pwd; %'/mnt/lab_48tb1/users/xiaohui/database/Sounding';

%Below is a list of variables containning information of stations.
% stationlist={'NKX'};
% lat_station=[32.85];
% lon_station=[-117.11];
% elev_station=[128];
% num_station=[72293];

stationlist={'YPPH'}; %{'NKX';'VBG';'OAK'};
lat_station=[-31.93]; %[32.85;34.75;37.73];
lon_station=[115.96]; %[-117.11;-120.56;-122.21];
elev_station=[20]; %[128;121;3.0];
num_station=[94610]; %[72293;72393;72493];

%Set time.
for yy=2014:2017
    for mm=1:12
        Year=num2str(yy); %{'2014';'2014'};
        Month=num2str(mm,'%02i'); %{'01';'02'};

        % Year={'2013';'2013';'2013';'2013';'2013';'2014';'2014';'2014';'2014';'2014';'2016'};
        % Month={'05';'06';'07';'08';'09';'05';'06';'07';'08';'09';'06'};
        % Start={'0100';'0100';'0100';'0100';'0100';'0100';'0100';'0100';'0100';'0100';'0100'};

        % Year={datestr(now,'yyyy')};
        % Month={datestr(now,'mm')};
        % Start={'0100'};
        % End= eomday(str2num(Year{1}),str2num(Month{1}));

        %Set header.
        %The HGHT in sounding data is geopotential height.
        header='PRES [hPa],HGHT [m],TEMP [C],DWPT [C],RELH [%],MIXR [g/kg],DRCT [deg],SKNT [knot],THTA [K],THTE [K],THTV [K], PW [mm], 1000 hPa to 500 hPa thickness [m]';
        %PRES in hPa; HGHT in m; THTA is Potential Temperature in K; MIXR is Mixing Ratio in g/kg; 

        %Firstly, loop through stationlist.
        for i=1:length(stationlist)

            %Secondly, loop through Year.
            for j=1:length(Year)

                End = eomday(yy,mm); %eomday(str2num(Year{j}),str2num(Month{j}));
                End = [num2str(End),'12'];

                %Get timerange.
                timerange=datenum(yy,mm,1,0,0,0):datenum(0,0,0,12,0,0):...
                          datenum(yy,mm,str2num(End(1:2)),str2num(End(3:4)),0,0);
%                 timerange=datenum(str2num(Year{j}),str2num(Month{j}),1,0,0,0):datenum(0,0,0,12,0,0):...
%                           datenum(str2num(Year{j}),str2num(Month{j}),str2num(End(1:2)),str2num(End(3:4)),0,0);

                %Get fileName and pathName to read.
                fileName=[];
                pathName=[];
                fileName=[num2str(num_station),'_',Year,'_',Month,'_','0100','_',End];
                pathName=fullfile(dirName,stationlist{i},'raw',fileName);

                %Read file using function read_mixed_csv.
                raw=[];
                raw=read_mixed_csv(pathName,' ');

                %Get the raw_temp is just the 1 column of all the information
                %contained in file specified by pathName.
                fid=fopen(pathName,'r');   
                raw_temp=[];     
                lineIndex=1;               
                nextLine=fgetl(fid);       
                while ~isequal(nextLine,-1)        
                    raw_temp{lineIndex,1}=nextLine;  
                    lineIndex=lineIndex+1;          
                    nextLine=fgetl(fid);           
                end
                fclose(fid); 
                %Remove empty cells, if needed.        
                raw_temp=raw_temp(1:lineIndex-1);  

                %Get row index of where the table head (PRES   HGHT   TEMP   DWPT
                %RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV) and tail
                %(</PRE><H3>Station information and sounding indices</H3><PRE>) is.
                index_head=[];
                index_head=find(strcmp(raw_temp,'   PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV')==1);
                %Add 3 to index_head to get the row index of the first row of data.
                index_head=index_head+3;
                index_tail=[];
                index_tail=find(strcmp(raw_temp,'</PRE><H3>Station information and sounding indices</H3><PRE>')==1);
                %Subtract 1 from index_tail to get the row index of the last row of
                %data.
                index_tail=index_tail-1;
                %Get row index of where the PW and 1000 hPa to 500 hPa thickness is.
                index_PW = strmatch('Precipitable water [mm] for entire sounding: ', raw_temp);
                index_thickness = strmatch('              1000 hPa to 500 hPa thickness: ', raw_temp);

                %Get fileName and pathName to write.
                fileName=[];
                pathName=[];
                fileName=[num2str(num_station(i)),'_',Year,'_',Month,'_','0100','_',End,'.csv'];
                pathName=fullfile(dirName,stationlist{i},'csv',fileName);

                %Write the Data into a file.
                fid=fopen(pathName,'w+');

                %Thirdly, loop through index.
                %Sometimes, either PW or 1000 hPa to 500 hPa thickness is missing
                %and cause code to crash, add counter and ensure 
                count_thickness = 0;
                count_PW = 0;
                for k=1:length(index_head)                        

                    temp=raw_temp{index_head(k)-6};
                    if strcmp(stationlist{i},'NKX')
                        time=datenum([temp(45:46),'-',temp(48:50),'-',temp(52:55)])+datenum(0,0,0,str2num(temp(41:42)),0,0);
                    elseif strcmp(stationlist{i},'VBG')
                        time=datenum([temp(50:51),'-',temp(53:55),'-',temp(57:60)])+datenum(0,0,0,str2num(temp(46:47)),0,0);
                    elseif strcmp(stationlist{i},'OAK')
                        time=datenum([temp(47:48),'-',temp(50:52),'-',temp(54:57)])+datenum(0,0,0,str2num(temp(43:44)),0,0);
                    elseif strcmp(stationlist{i},'YPPH')
                        time=datenum([temp(50:51),'-',temp(53:55),'-',temp(57:60)])+datenum(0,0,0,str2num(temp(46:47)),0,0);
                    end

                    %Write the date.            
                    fprintf(fid,'%s\n',datestr(time,'yyyymmdd_HH'));
                    %Write header.
                    fprintf(fid,'%s\n',header);

                    %Initialize Data to be written.
                    Data=NaN(index_tail(k)-index_head(k)+1,11);

                    %Read data.
                    temp=[];
                    temp=raw(index_head(k):index_tail(k),:);

                    %Loop through rows of temp.
                    for row=1:length(temp(:,1))

                        %The count is used as a counter.
                        count=0;
                        %Loop through columns of temp.
                        for col=1:length(temp(1,:))

                            %Determine whether temp{row,col} is empty or not, this
                            %method is not prefect, as sometimes, water vapor
                            %mixing ratio and RH data is missing after QVAPOR
                            %becomes 0. 
                            if isempty(temp{row,col})~=1

                                count=count+1;
                                Data(row,count)=str2num(temp{row,col});

                            end

                        end

                        %Write Data.
                        if row == 1
                            fprintf(fid,'%1f,%d,%.1f,%.1f,%d,%.2f,%d,%d,%.1f,%.1f,%.1f',Data(row,1),Data(row,2),Data(row,3),Data(row,4),...
                                    Data(row,5),Data(row,6),Data(row,7),Data(row,8),Data(row,9),Data(row,10),Data(row,11));                        

                            %Write Precipitable water (PW [mm]),  1000 hPa to 500 hPa thickness [m]  
                            PW = nan;
                            thickness_500mb = nan;
                            if k < length(index_head)
                                if count_PW <= length(index_PW) & index_PW(count_PW+1) > index_tail(k) & index_PW(count_PW+1) < index_head(k+1)
                                    count_PW = count_PW +1;
                                    PW = str2num(raw_temp{index_PW(count_PW)}(46:end));
                                end
                                if index_thickness(count_thickness+1) > index_tail(k) && index_thickness(count_thickness+1) < index_head(k+1)
                                    count_thickness = count_thickness +1;
                                    thickness_500mb = str2num(raw_temp{index_thickness(count_thickness)}(46:end));
                                end

                            else
                                if count_PW <= length(index_PW) & index_PW(count_PW+1) > index_tail(k) 
                                    count_PW = count_PW +1;
                                    PW = str2num(raw_temp{index_PW(count_PW)}(46:end));
                                end
                                if index_thickness(count_thickness+1) > index_tail(k)
                                    count_thickness = count_thickness +1;
                                    thickness_500mb = str2num(raw_temp{index_thickness(count_thickness)}(46:end));
                                end

                            end

                            %Check whether the values of PW or thickness_500mb is
                            %normal.
                            if PW < 1
                               PW = nan;                       
                            end
                            if thickness_500mb < 2000
                                thickness_500mb = nan;
                            end

                            fprintf(fid,',%.1f,%d\n', PW, thickness_500mb);

                        else
                            fprintf(fid,'%1f,%d,%.1f,%.1f,%d,%.2f,%d,%d,%.1f,%.1f,%.1f\n',Data(row,1),Data(row,2),Data(row,3),Data(row,4),...
                                    Data(row,5),Data(row,6),Data(row,7),Data(row,8),Data(row,9),Data(row,10),Data(row,11));                        
                        end

                    end            

                end                         

                fclose(fid);
                fprintf('%s is written!\n',pathName);

            end

        end
        
    end %mm
end %yy
