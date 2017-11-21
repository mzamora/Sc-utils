function [DT_max,numofInv,varargout]=TMP_Inversion_Strength_Cal_V1(T,hght,hght_sfc)
% [DT_max,numofInv,varargout]=TMP_Inversion_Strength_Cal_V1(T,hght,hght_sfc)
% The input includes:
% T: temperature (unit: K) of 1 dimension (z)
% hght: height (unit: km) of 1 dimension (z)
% hght_sfc: terrain height (unit: m) of 1 dimension (z)
% Written by Xiaohui Zhong, SRAF UCSD http://solar.ucsd.edu
% Modified my Monica Zamora

%Convert ground elevation from m to km.
hght_sfc=hght_sfc./1000;

%The p_Crit (hPa) and H_Crit (km) specifies the pressure and altitude to
%which from the surface the temperature is examined.
p_Crit=700;
H_Crit=3;
% H_Crit=2;

%Initialize outputs.
DT_max=[];
DT_DZ_max=[];
numofInv=[];
%The hght_top and hght_base is the height of top and base of the
%temperature inversion layer.
%The eta_top and eta_base and  is the eta level of top and base of the
%temperature inversion layer.
hght_top=[];
hght_base=[];
eta_top=[];
eta_base=[];    

%Initialized.
%The eta_top and eta_end are defined as the eta levels for each of
%temperature inversion layers within 3 km.
eta_top_temp=[];
eta_base_temp=[];
DT=[];
DT_DZ=[];
T_PROF=[];
HGHT=[];
if size(T,2)>size(T,1)
    T_PROF=T';
    HGHT=hght';
else
    T_PROF=T;
    HGHT=hght;
end
%Restrict height within H_Crit km.
index=[];
index=find(HGHT<H_Crit);
HGHT=HGHT(index);
T_PROF=T_PROF(index);

%Calculate whether temperature is increasing or decreasing with
%height.
dT=sign(T_PROF(2:end)-T_PROF(1:end-1));
%         dT(end+1)=dT(end);

%Find first positive dT.
INDICES=[];
INDICES=find(dT>0);

if isempty(INDICES)~=1

    % find longest string of consecutive indices
    dI=[];
    dI=INDICES(2:end)-INDICES(1:end-1);

            
%     %{
    %%%For the case when a small decrease in temperature with
    %%%height within a inversion causes two inversions, it is
    %%%necessary to treat them as one inverison when the lower inversion is not surface inversion.     
    %Determine whether there is a surface inversion to avoid treating surface
    %and non-surface inversion as one inversion.
    if HGHT(INDICES(1))-hght_sfc<=0.05
        start=2;
    else
        start=1;        
    end
    %The threshold is used to determine whether the difference
    %is within the reasonable range.   
    threshold=0.5;
    temp=find(dI==2);            
    if isempty(temp)~=1        
        for m=start:length(temp)
            if abs(T_PROF(INDICES(temp(m))+1)-T_PROF(INDICES(temp(m)+1)))<=threshold
                INDICES(temp(m)+2:end+1)=INDICES(temp(m)+1:end);
                INDICES(temp(m)+1)=INDICES(temp(m))+1;

                temp=temp+1;                        

            end
        end
    end        

    %There are also cases when 2 grid points separating two
    %inversions, which should be treated as one.
    dI=INDICES(2:end)-INDICES(1:end-1);
    temp=find(dI==3);
    if isempty(temp)~=1
        for m=start:length(temp)
            if abs(T_PROF(INDICES(temp(m))+1)-T_PROF(INDICES(temp(m)+1)))<=threshold
                INDICES(temp(m)+3:end+2)=INDICES(temp(m)+1:end);
                INDICES(temp(m)+1)=INDICES(temp(m))+1;
                INDICES(temp(m)+2)=INDICES(temp(m))+2;

                temp=temp+2;                        

            end
        end
    end            
    dI=INDICES(2:end)-INDICES(1:end-1);            

    %There are also cases when 3 grid points separating two
    %inversions, which should be treated as one.
    dI=INDICES(2:end)-INDICES(1:end-1);
    temp=find(dI==4);
    if isempty(temp)~=1
        for m=start:length(temp)
            if abs(T_PROF(INDICES(temp(m))+1)-T_PROF(INDICES(temp(m)+1)))<=threshold
                INDICES(temp(m)+4:end+3)=INDICES(temp(m)+1:end);
                INDICES(temp(m)+1)=INDICES(temp(m))+1;
                INDICES(temp(m)+2)=INDICES(temp(m))+2;
                INDICES(temp(m)+3)=INDICES(temp(m))+3;

                temp=temp+3;                        

            end
        end
    end            
    dI=INDICES(2:end)-INDICES(1:end-1);
    %%%
    %}                
    
    SEG_INDEX=find(dI>1);
    N_SEGS=length(SEG_INDEX)+1;

    %Add the index of last temperature in the temperature inversion.
    if N_SEGS==1
        INDICES=[INDICES;INDICES(end)+1];
    elseif N_SEGS==2
        INDICES=[INDICES(1:SEG_INDEX);INDICES(SEG_INDEX)+1;INDICES(SEG_INDEX+1:end);INDICES(end)+1];
        SEG_INDEX=SEG_INDEX+1;
        
    elseif N_SEGS>2

        INDICES_temp=[INDICES(1:SEG_INDEX(1));INDICES(SEG_INDEX(1))+1];

        for n=1:N_SEGS-2
            INDICES_temp=[INDICES_temp;INDICES(SEG_INDEX(n)+1:SEG_INDEX(n+1));INDICES(SEG_INDEX(n+1))+1];
            SEG_INDEX(n)=SEG_INDEX(n)+n;
        end               

        INDICES_temp=[INDICES_temp;INDICES(SEG_INDEX(N_SEGS-1)+1:end);INDICES(end)+1];
        INDICES=INDICES_temp;
        SEG_INDEX(N_SEGS-1)=SEG_INDEX(N_SEGS-1)+n+1;

    end

    % find length of each segment SEG_LENGTHS.

    %Calculate inversion strength within each of inversion layers.
    if N_SEGS==1
        SEG_LENGTHS=length(INDICES);

        DT=T_PROF(INDICES(end))-T_PROF(INDICES(1));
        DT_DZ=(T_PROF(INDICES(end))-T_PROF(INDICES(1)))/(HGHT(INDICES(end))-HGHT(INDICES(1)));

        eta_top_temp=INDICES(end);
        eta_base_temp=INDICES(1);
    
    elseif N_SEGS==2
        SEG_LENGTHS(1)=length(INDICES(1:SEG_INDEX(1)));            
        SEG_LENGTHS(2)=length(INDICES(SEG_INDEX(1)+1:end));

        DT(1)=T_PROF(INDICES(SEG_INDEX(1)))-T_PROF(INDICES(1));
        DT(2)=T_PROF(INDICES(end))-T_PROF(INDICES(SEG_INDEX(1)+1));
        DT_DZ(1)=(T_PROF(INDICES(SEG_INDEX(1)))-T_PROF(INDICES(1)))/(HGHT(INDICES(SEG_INDEX(1)))-HGHT(INDICES(1)));
        DT_DZ(2)=(T_PROF(INDICES(end))-T_PROF(INDICES(SEG_INDEX(1)+1)))/(HGHT(INDICES(end))-HGHT(INDICES(SEG_INDEX(1)+1)));

        eta_top_temp(1)=INDICES(SEG_INDEX(1));
        eta_top_temp(2)=INDICES(end);
        eta_base_temp(1)=INDICES(1);
        eta_base_temp(2)=INDICES(SEG_INDEX(1)+1);
        
    elseif N_SEGS>2
        SEG_LENGTHS(1)=length(INDICES(1:SEG_INDEX(1)));
        
        DT(1)=T_PROF(INDICES(SEG_INDEX(1)))-T_PROF(INDICES(1));
        DT_DZ(1)=(T_PROF(INDICES(SEG_INDEX(1)))-T_PROF(INDICES(1)))/(HGHT(INDICES(SEG_INDEX(1)))-HGHT(INDICES(1)));

        eta_top_temp(1)=INDICES(SEG_INDEX(1));
        eta_base_temp(1)=INDICES(1);
        
        for n=2:N_SEGS-1
            SEG_LENGTHS(n)=length(INDICES(SEG_INDEX(n-1)+1:SEG_INDEX(n)));
            
            DT(n)=T_PROF(INDICES(SEG_INDEX(n)))-T_PROF(INDICES(SEG_INDEX(n-1)+1));
            DT_DZ(n)=(T_PROF(INDICES(SEG_INDEX(n)))-T_PROF(INDICES(SEG_INDEX(n-1)+1)))/(HGHT(INDICES(SEG_INDEX(n)))-HGHT(INDICES(SEG_INDEX(n-1)+1)));
            
            eta_top_temp=[eta_top_temp,INDICES(SEG_INDEX(n))];
            eta_base_temp=[eta_base_temp,INDICES(SEG_INDEX(n-1)+1)];
            
        end

        SEG_LENGTHS(N_SEGS)=length(INDICES(SEG_INDEX(N_SEGS-1)+1:end));
        
        DT(N_SEGS)=T_PROF(INDICES(end))-T_PROF(INDICES(SEG_INDEX(N_SEGS-1)+1));
        DT_DZ(N_SEGS)=(T_PROF(INDICES(end))-T_PROF(INDICES(SEG_INDEX(N_SEGS-1)+1)))/(HGHT(INDICES(end))-HGHT(INDICES(SEG_INDEX(N_SEGS-1)+1)));
        
        eta_top_temp=[eta_top_temp,INDICES(end)];
        eta_base_temp=[eta_base_temp,INDICES(SEG_INDEX(N_SEGS-1)+1)];

    end

    %Get the inversion strength of the strongest temperature
    %inversion.
    DT_max=max(DT);
    DT_DZ_max=max(DT_DZ);
    numofInv=N_SEGS;     

    if numofInv==1
        eta_top=eta_top_temp;
        eta_base=eta_base_temp;
        hght_top=hght(eta_top_temp);
        hght_base=HGHT(eta_base_temp); 

    elseif numofInv>1
        %Find the strongest temperature inversion except for the surface
        %inversion.
        index=find(DT==DT_max);

        if index(1)==1 && HGHT(eta_base_temp(1))-hght_sfc<=0.05 && HGHT(eta_top_temp(1))-hght_sfc<=0.1                    
           %Surface (radiation) inversion is strongest and there is
           %more than surface inversion.

           DT(1)=[];
           DT_DZ(1)=[];
           eta_base_temp(1)=[];
           eta_top_temp(1)=[];

           DT_max=max(DT);
           DT_DZ_max=max(DT_DZ);
           index=find(DT==DT_max);

        end
              
        eta_top=eta_top_temp(index(1));
        eta_base=eta_base_temp(index(1));
        hght_top=hght(eta_top_temp(index(1)));                
        hght_base=HGHT(eta_base_temp(index(1)));        
        
    end
    
    if DT_max<3
        %The inversion is not strong enough to be called inversion
        %here.
        eta_top=NaN;
        eta_base=NaN;
        hght_top=NaN;
        hght_base=NaN;
    end        

else
    eta_top=NaN;
    eta_base=NaN;
    hght_top=NaN;
    hght_base=NaN;
    DT_max=0;
    numofInv=0;

end


%% Checking inversion top height with DTDz<5
z=hght;  
dTdz=diff(T)./diff(z);
istar=eta_base+find(dTdz(eta_base+1:eta_top+1)<5,1,'first');
hght_top2=z(istar);
if hght_top2<hght_top
    hght_top=hght_top2;
    eta_top=istar;
    DT_max=T(eta_top)-T(eta_base);
    if DT_max<3
        %The inversion is not strong enough to be called inversion
        %here.
        eta_top=NaN;
        eta_base=NaN;
        hght_top=NaN;
        hght_base=NaN;
    end      
end  

%% Output variables
if nargout>=6
    varargout{1}=hght_top;
    varargout{2}=hght_base;    
    varargout{3}=eta_top;
    varargout{4}=eta_base;
end

end
