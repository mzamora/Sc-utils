function  def=def_cnstnts
%  Description: Set of atmospheric constants, surface sensible, latent heat  
%               flux, radiation parameter values
%
%  Author: Mohamed Ghonima
%% Specifiy constants
def.R=287.05; %Ra
def.Rm=461.5; % Rv
def.ep=def.R/def.Rm;
def.cp=1005;
def.rcp=def.R/def.cp;
def.cpr=def.cp/def.R;
def.p00=1.0*10^5;
def.alvl=2.5*10^6;
def.g=9.81;
def.psrf=1.0177e05;
def.sigma=5.67*10^-8;
def.th00=289;
