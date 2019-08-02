function varargout=Get_sounding_Var(time,Time_End,pathName,ListofVar)
% [p z tp td RH w winddir windspeed theta theta_e theta_v]=Get_sounding_Var(time,Time_End,pathName,ListofVar)
% Developed by Xiaohui Zhong, UCSD SRAF http://solar.ucsd.edu
% It reads csv files downloaded from UWyoming's sounding database
% Slightly modified by Monica Zamora, UCSD SRAF (Nov 2017) http://solar.ucsd.edu


% addpath(genpath('/home/xiaohui/WRF/m_programe_Linux'));
% addpath(genpath('/home/xiaohui/REST2/m_programe'));
VarList={'PRES';'HGHT';'TEMP';'DWPT';'RELH';'MIXR';'WDIR';'WSPD';'THTA';'THTE';'THTV'};

%Read file using function read_mixed_csv.
raw=[];
raw=read_mixed_csv(pathName,',');

%Get row index where the data corresponding to time is in Sounding data
%file.
index_Start=[];
index_Start=find(strcmp(datestr(time,'yyyymmdd_HH'),raw)==1)+2;
index_End=[];
if abs(time-Time_End)<datenum(0,0,0,0,1,0)
    index_End=length(raw(:,1));
else
    index_End=find(strcmp(datestr(time+datenum(0,0,0,12,0,0),'yyyymmdd_HH'),raw)==1)-1;
    
    if isempty(index_End)==1
        %This happens because Sounding measurement is only taken at 12 UTC
        %on each day.
        index_End=find(strcmp(datestr(time+datenum(0,0,0,24,0,0),'yyyymmdd_HH'),raw)==1)-1;        
    end
    
end

%Loop through the list of variables given in the input.
for i=1:length(ListofVar)
    Var=[];
    index=find(strcmp(VarList,ListofVar{i})==1);
    if isempty(index)~=1
        Var=str2num(sprintf('%s\n',raw{index_Start:index_End,index}));
        switch ListofVar{i}            
            case {'TEMP';'DWPT'}
                Var=Var+273.15;
                
            case 'WSPD'
                %Convert wind speed from knot to m/s.
                Var=Var*0.514;
                
        end
    else
        return;
        disp('Wrong input variable name!');
    end
    
    varargout{i}=Var;
end                 

%% Check that there aren't duplicated values in the BL
p=varargout{1};
tp=varargout{3};
z=varargout{2};
varargout{12}=z(p==1000); %geopotential height
varargout{13}=z(p==500); %geopotential height

idelete=[];
for istar=1:length(tp)-1 %go through all levels
    if (tp(istar)==tp(istar+1)) %the temp is the same up or down
         idelete=[idelete, istar']; %add this index to the ones to be deleted
    end
end

for i=1:length(VarList) %go through the indices we want to delete
    varargout{i}(idelete)=[];
end

%% Sanity check: Remove all rows with nan values
for j=1:length(ListofVar) %go through all variables
    f=isnan(varargout{j}); %find nans 
    for i=1:length(ListofVar)
        varargout{i}(f)=[]; %remove nan rows everywhere
    end
end

end
