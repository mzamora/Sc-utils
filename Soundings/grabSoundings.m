%grab soundings from UWyoming website

for yy=2016
    for mm=1:12
        stationnumber='85442';
        station='SCFA';
        End=[num2str(eomday(yy,mm)),'12'];
        out=[station,'/raw/',stationnumber,'_',num2str(yy),'_',num2str(mm,'%02i'),'_0100_',End];
        site=['http://weather.uwyo.edu/cgi-bin/sounding?region=samer&TYPE=TEXT%3ALIST&YEAR=', ...
            num2str(yy),'&MONTH=',num2str(mm,'%02i'),'&FROM=0100&TO=',End,'&STNM=',stationnumber];
        options = weboptions('Timeout',Inf);
        websave(out,site,options);
    end
end