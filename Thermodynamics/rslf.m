function [rslf]=rslf(p,t)

% This function calculates the water saturation vapor mixing ratio as a
% function of temperature and pressure

tmelt  = 273.16; % Triple point temperature

% A six order poynomial with the folloqing cosntants is  used to calculate
% saturation pressure as a function of Temp
c0=0.6105851e+03; c1=0.4440316e+02;c2=0.1430341e+01; c3=0.2641412e-01;
c4=0.2995057e-03; c5=0.2031998e-05;c6=0.6936113e-08; c7=0.2564861e-11;
c8=-.3704404e-13;

x=max(-80,t-tmelt);
esl=c0+x.*(c1+x.*(c2+x.*(c3+x.*(c4+x.*(c5+x.*(c6+x.*(c7+x.*c8)))))));

% Calculate saturation mixing ratio
rslf=.622*esl./(p-esl);