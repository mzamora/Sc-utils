%% get_raw_sounding.m for a station, and years wanted

%% Stations
%Perth: YPPH 94610
stname='YPPH'; stnumber='94610'; region='pac';
%Learmonth: YPLM 94302
%stname='YPLM'; stnumber='94302'; region='pac';
% 85442 SCFA Antofagasta
%stnumber='85442'; stname='SCFA'; region='samer';
%Pointe-Noire 64400 FCPP
%stname='FCPP'; stnumber='64400'; region='africa';

%% Obtain html
url0=['http://weather.uwyo.edu/cgi-bin/sounding?region=',region,'&TYPE=TEXT%3ALIST&'];
for year=2014:2017
    for month=1:12
        dend=[num2str(eomday(year,month)),'12'];
        url=[url0,'YEAR=',num2str(year),'&MONTH=',num2str(month,'%02i'),'&FROM=0100&TO=',dend,'&STNM=94610'];
        filename=[stname,'/raw/',stnumber,'_',num2str(year),'_',num2str(month,'%02i'),'_0100_',dend];
        try
            websave(filename,url);
        catch
            fprintf(['Could not get web file',num2str(year),'_',num2str(month,'%02i'),' \n'])
        end
    end
end
        