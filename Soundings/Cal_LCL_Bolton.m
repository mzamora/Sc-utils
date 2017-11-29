function [T_LCL2,P_LCL2,H_LCL2]=Cal_LCL_Bolton(TSK,PSFC,QVG,QSG,HGT,Tdew)
% [T_LCL2,P_LCL2,H_LCL2]=Cal_LCL_Bolton(TSK,PSFC,QVG,QSG,HGT,Tdew)
% There exist several different formulas for calculating the temperature and
% pressure of the LCL.
% This function is based on the 1980_ Bolton_ The Computation of Equivalent
% Potential Temperature.
% The units of TSK, PSFC and QVG must be K, Pa and kg/kg.
% Written by Xiaohui Zhong, SRAF UCSD http://solar.ucsd.edu

%Define some constants
%The Rd is specific gas constant, dry air (J/(kg-K)).
Rd=287.04; 
%The Rv is specific gas constant (J/(kg-K)).
Rv=461.5; 
%The L is lapse rate (K/m).
L=-6.5*10^(-3);
%The g is gravitional constant (m/s^2).
g=9.80665;

%Get TSK_temp in C degrees.
TSK_temp=TSK-273.15;

%Initialize.
%D_SFCn are calculated using different methods or equations.
D_SFC1=[];
D_SFC2=[];
D_SFC3=[];
D_SFC4=[];

%T_LCLn, P_LCLn and H_LCLn are calculated using different D_SFCn.
T_LCL1=[];
T_LCL2=[];
P_LCL1=[];
P_LCL2=[];
H_LCL1=[];
H_LCL2=[];

w_s=[];

for i=1:size(PSFC,1)
   for j=1:size(PSFC,2) 
        %Calculate saturation water vapor pressure using the clausius
        %clapeyron equation.
        %The returned values of e_s are in hPa.
        e_s=[];
        tc=TSK(i,j)-273.15;
        e_s=610.94*exp(17.652*tc/(tc+243.04));
       
        %Calculate P_d (dry partial pressure) in Pa.
        PSFC_d=[];
        PSFC_d=PSFC(i,j)-e_s;

        %Calculate w_s (saturation water vapor mixing ratio) in kg/kg.        
        w_s(i,j)=Rd/Rv.*(e_s./PSFC_d);
        
        %Calculate relative humidity RH_SFC in %.
        RH_SFC=[];
        if isempty(QSG)==1
            RH_SFC=100*QVG(i,j)./w_s(i,j);        
        else
            RH_SFC=100*QVG(i,j)./QSG(i,j);        
        end
        
        %%
        %Calculate dewpoint temperature using different formulas of
        %functions.
        
        %%%
        %Based on dewpoint formulas from
        %http://ag.arizona.edu/azmet/dewpoint.html.

        %Calculate B, an intermediate value (no units).
        B=[];
        B=(log(RH_SFC/100)+((17.27*TSK_temp(i,j))/(237.3+TSK_temp(i,j))))/17.27;
        
        %Calculate dew point temperature in C degrees.
%         D_SFC1_temp(i,j)=(237.3*B)./(1-B);
%         %Convert D_SFC1 from C degrees to K.
%         D_SFC1(i,j)=D_SFC1_temp(i,j)+273.15;                
%         
%         %%%
%         %Calculate dew point temperature using the function
%         %convert_humidity and 'Bolton1980' as the method.
%         D_SFC2(i,j)=convert_humidity(PSFC(i,j),TSK(i,j),QVG(i,j),'mixing ratio','dew point temperature','Bolton1980');
%         
%         %%%
%         %Based on dewpoint formulas from 
%         %http://andrew.rsmas.miami.edu/bmcnoldy/humidity_conversions.pdf. 
%         b=243.04;
%         a=17.625;
%         D_SFC3(i,j)=b*(log(RH_SFC/100)+a*TSK_temp(i,j)/(b+TSK_temp(i,j)))/(a-log(RH_SFC/100)-a*TSK_temp(i,j)/(b+TSK_temp(i,j)));
%         D_SFC3(i,j)=D_SFC3(i,j)+273.15;
%         
%         %%%
%         %Based on a simpler calculation that gives an approximation of dew
%         %point temperature from
%         %http://iridl.ldeo.columbia.edu/dochelp/QA/Basic/dewpoint.html.
%         D_SFC4(i,j)=TSK(i,j)-((100 - RH_SFC)/5);
        %%
        %Calculate temperature at lifting condensation level (LCL).
%         T_LCL1(i,j)=1/(1/(D_SFC1(i,j)-56)+log(TSK(i,j)/D_SFC1(i,j))/800)+56;
        T_LCL2(i,j)=1/(1/(Tdew(i,j)-56)+log(TSK(i,j)/Tdew(i,j))/800)+56;
        
        %Calculate pressure at LCL based on the definition of potential
        %temperature (Poisson's equation).
%         P_LCL1(i,j)=PSFC(i,j)*(T_LCL1(i,j)/TSK(i,j))^3.5;
        P_LCL2(i,j)=PSFC(i,j)*(T_LCL2(i,j)/TSK(i,j))^3.5;
                        
%         %Calculate height at LCL based on the hyposmetric formula from
%         %http://keisan.casio.com/exec/system/1224585971.
%         H_LCL1(i,j)=((PSFC(i,j)/P_LCL1(i,j))^(1/5.257)-1)*T_LCL1(i,j)/0.0065+HGT(i,j);
%         H_LCL2(i,j)=((PSFC(i,j)/P_LCL2(i,j))^(1/5.257)-1)*T_LCL2(i,j)/0.0065+HGT(i,j);
        
        %Calculate height at LCL based on the hyposmetric equation from
        %http://en.wikipedia.org/wiki/Hypsometric_equation.
        %The mean temperature between the surface and LCL are calculated as
        %mean value between TSK and T_LCL.
%         H_LCL1(i,j)=Rd*(TSK(i,j)+T_LCL1(i,j))/2/g*log(PSFC(i,j)/P_LCL1(i,j))+HGT(i,j);
        H_LCL2(i,j)=Rd*(TSK(i,j)+T_LCL2(i,j))/2/g*log(PSFC(i,j)/P_LCL2(i,j))+HGT(i,j);
        
%         %Calculate height at LCL based on the equation from
%         %http://psas.pdx.edu/RocketScience/PressureAltitude_Derived.pdf.
%         H_LCL1(i,j)=TSK(i,j)/L*((P_LCL1(i,j)/PSFC(i,j))^(-L*Rd/g)-1)+HGT(i,j);
%         H_LCL2(i,j)=TSK(i,j)/L*((P_LCL2(i,j)/PSFC(i,j))^(-L*Rd/g)-1)+HGT(i,j);
        
   end
    
end


end
