function [istar,jstar]=find_NAM_indices(lat,lon)
% find the indices for the latitudes and longitudes given
% needs the NAM csv

load NAM_lon.csv
load NAM_lat.csv

n=length(lat);

for i=1:n
    deltamin=100;

    %Find the indices that minimizes the euclidean distance to the desired coordinates
    for ii=1:size(NAM_lon,1)
        for jj=1:size(NAM_lon,2)
            delta=abs(lat(i)-NAM_lat(ii,jj))^2+abs(lon(i)-NAM_lon(ii,jj))^2;
            
            if delta<deltamin
                deltamin=delta;
                istar(i)=ii; %y
                jstar(i)=jj; %x
            end
        end
    end
end
