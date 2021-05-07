function Tw=T_wetbulb(T,RH)
% Tw=T_wetbulb(T,RH)
% all temp in Celsius, all RH in [0,100]

% wet bulb. stull, 2011 (doi 10.1175/JAMC-D-11-0143.1)
Tw=T.*atan(0.151977 * sqrt( RH + 8.313659 )) + atan(T + RH) - atan(RH-1.676331) + ...
  0.00391838*(RH.^1.5)*atan(0.023101*RH)-4.686035;

end %function
