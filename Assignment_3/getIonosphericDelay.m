function [ionosphericDelay] = getIonosphericDelay(TECU,Zenit_angle)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% TECU=5.3;
%constant parameter
Re=6371000; %m
f=1575.42*1e6;
hi=350000; %m

% start function
TEC=TECU*10^16;

OF=(1-((Re*sind(Zenit_angle)/(Re+hi)).^2)).^(-0.5)
ionosphericDelay=(40.3/(f^2))*TEC*OF

end