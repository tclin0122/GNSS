function [ionosphericDelay] = getIonosphericDelay(TECU,Zenit_angle,Height)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% TECU=5.3;
%constant parameter
Re=6371000; %m
f=1575.42*1e6;
% start function
TEC=TECU*10^16;
zen=90-el;
OF=(1-((Re*sind(Zenit_angle)/(Re+Height)).^2)).^(-0.5)
ionosphericDelay=(40.3/(f^2))*TEC*OF
end