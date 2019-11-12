function Tdp=T_dewpoint(T,RH)
% all temp in Celsius, all RH in [0,100]

% dew point. bolton, 1980
b=18.678; c=257.14; %deg C
gamma=log(RH/100)+(b*T)./(c*T);
Tdp=c*gamma./(b-gamma);

end %function
