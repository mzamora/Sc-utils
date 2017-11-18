function [T_LCL2,P_LCL2,H_LCL2,varargout]=Cal_LCL_Barnes(TSK,PSFC,QVG,QSG,HGT,Tdew)
% [T_LCL2,P_LCL2,H_LCL2,varargout]=Cal_LCL_Barnes(TSK,PSFC,QVG,QSG,HGT,Tdew)
% There exist several different formulas for calculating the temperature and
% pressure of the LCL.
% This function is based on the 1968_ Stanley L. Barnes_ An Empirical
% Shortcut to the Calculation of Temperature and Pressure at the Lifted
% Condensation Level.
% The units of TSK, PSFC, QVG, QSG and HGT must be K, Pa, kg/kg, kg/kg and
% m. Written by Xiaohui Zhong, SRAF UCSD http://solar.ucsd.edu

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

for m=1:size(PSFC,3)
    
    for i=1:size(PSFC,1)

       for j=1:size(PSFC,2) 
            %Calculate saturation water vapor pressure using the clausius
            %clapeyron equation.
            %The returned values of e_s are in hPa.
            e_s=[];
            e_s=clausius_clapeyron(TSK(i,j,m));
            %Convert e_s from hPa to Pa.
            e_s=e_s*100;

            %Calculate P_d (dry partial pressure) in Pa.
            PSFC_d=[];
            PSFC_d=PSFC(i,j,m)-e_s;

            %Calculate w_s (saturation water vapor mixing ratio) in kg/kg.        
            w_s(i,j,m)=Rd/Rv.*(e_s./PSFC_d);

            %Calculate relative humidity RH_SFC in %.
            RH_SFC=[];
            if isempty(QSG)==1
                RH_SFC=100*QVG(i,j,m)./w_s(i,j,m);        
            else
                RH_SFC=100*QVG(i,j,m)./QSG(i,j,m);        
            end

            %%
            %Calculate dewpoint temperature using different formulas of
            %functions.

            %%%
            %Based on dewpoint formulas from
            %http://ag.arizona.edu/azmet/dewpoint.html .

            %Calculate B, an intermediate value (no units).
            B=[];
            B=(log(RH_SFC/100)+((17.27*TSK_temp(i,j,m))/(237.3+TSK_temp(i,j,m))))/17.27;

%             %Calculate dew point temperature in C degrees.
%             D_SFC1(i,j,m)=(237.3*B)./(1-B);
%             %Convert D_SFC1 from C degrees to K.
%             D_SFC1(i,j,m)=D_SFC1(i,j,m)+273.15;                
% 
%             %%%
%             %Calculate dew point temperature using the function
%             %convert_humidity and 'Bolton1980' as the method.
%             D_SFC2(i,j,m)=convert_humidity(PSFC(i,j,m),TSK(i,j,m),QVG(i,j,m),'mixing ratio','dew point temperature','Bolton1980');
% 
%             %%%
%             %Based on dewpoint formulas from 
%             %http://andrew.rsmas.miami.edu/bmcnoldy/humidity_conversions.pdf . 
%             b=243.04;
%             a=17.625;
%             D_SFC3(i,j,m)=b*(log(RH_SFC/100)+a*TSK_temp(i,j,m)/(b+TSK_temp(i,j,m)))/(a-log(RH_SFC/100)-a*TSK_temp(i,j,m)/(b+TSK_temp(i,j,m)));
%             D_SFC3(i,j,m)=D_SFC3(i,j,m)+273.15;
% 
%             %%%
%             %Based on a simpler calculation that gives an approximation of dew
%             %point temperature from
%             %http://iridl.ldeo.columbia.edu/dochelp/QA/Basic/dewpoint.html .
%             D_SFC4(i,j,m)=TSK(i,j,m)-((100 - RH_SFC)/5);
            %%
            %Calculate temperature at lifting condensation level (LCL).
%             T_LCL1(i,j,m)=D_SFC1(i,j,m)-(0.001296*D_SFC1(i,j,m)+0.1963)*(TSK(i,j,m)-D_SFC1(i,j,m));
%             T_LCL2(i,j,m)=D_SFC2(i,j,m)-(0.001296*D_SFC2(i,j,m)+0.1963)*(TSK(i,j,m)-D_SFC2(i,j,m));
            T_LCL2(i,j,m)=Tdew(i,j,m)-(0.001296*Tdew(i,j,m)+0.1963)*(TSK(i,j,m)-Tdew(i,j,m));

            %Calculate pressure at LCL based on the definition of potential
            %temperature (Poisson's equation).
%             P_LCL1(i,j,m)=PSFC(i,j,m)*(T_LCL1(i,j,m)/TSK(i,j,m))^3.5;
            P_LCL2(i,j,m)=PSFC(i,j,m)*(T_LCL2(i,j,m)/TSK(i,j,m))^3.5;

    %         %Calculate height at LCL based on the hyposmetric formula from
    %         %http://keisan.casio.com/exec/system/1224585971.
    %         H_LCL1(i,j,m)=((PSFC(i,j,m)/P_LCL1(i,j,m))^(1/5.257)-1)*T_LCL1(i,j,m)/0.0065+HGT(i,j,1);
    %         H_LCL2(i,j,m)=((PSFC(i,j,m)/P_LCL2(i,j,m))^(1/5.257)-1)*T_LCL2(i,j,m)/0.0065+HGT(i,j,1);

            %Calculate height at LCL based on the hyposmetric equation from
            %http://en.wikipedia.org/wiki/Hypsometric_equation .
            %The mean temperature between the surface and LCL are calculated as
            %mean value between TSK and T_LCL.
%             H_LCL1(i,j,m)=Rd*(TSK(i,j,m)+T_LCL1(i,j,m))/2/g*log(PSFC(i,j,m)/P_LCL1(i,j,m))+HGT(i,j,1);
            H_LCL2(i,j,m)=Rd*(TSK(i,j,m)+T_LCL2(i,j,m))/2/g*log(PSFC(i,j,m)/P_LCL2(i,j,m))+HGT(i,j,1);

    %         %Calculate height at LCL based on the equation from
    %         %http://psas.pdx.edu/RocketScience/PressureAltitude_Derived.pdf.
    %         H_LCL1(i,j,m)=TSK(i,j,m)/L*((P_LCL1(i,j,m)/PSFC(i,j,m))^(-L*Rd/g)-1)+HGT(i,j,1);
    %         H_LCL2(i,j,m)=TSK(i,j,m)/L*((P_LCL2(i,j,m)/PSFC(i,j,m))^(-L*Rd/g)-1)+HGT(i,j,1);

       end

    end

end

if nargout==4
    varargout{1}=D_SFC2;
end

end
