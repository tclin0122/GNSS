function [delaySaas, delayCorr] = getTroposfericDelay(Zenit_angle, Height)

% Parameters
P = 1013;           % [mBar,hPa] 
T = 273.15;     % [Kelvin]
RH = 0.5;

% Read data files
B_data = readtable('B_data.csv');
d_R_data = readmatrix('d_R_data.csv');

% Simulation parameters
interp_resolution = 0.05; 

% B data
height = B_data.Height;
B = B_data.B;
old_height = d_R_data(1,2:end);
ipl_height = min(old_height):interp_resolution:max(old_height);
B_interp = interp1(height, B, ipl_height);
B_table = array2table([transpose(ipl_height) transpose(B_interp)]);
B_table.Properties.VariableNames(1:2) = {'Height','B'};

% d_R data
degrees = d_R_data(2:end,1);
values = d_R_data(2:end,2:end);
[X,Y] = meshgrid(old_height,degrees);
degrees_interp = min(degrees):0.5:max(degrees);
[Xq,Yq] = meshgrid(ipl_height,degrees_interp);
d_R_interp = interp2(X,Y,values,Xq,Yq);
d_R_table = array2table([[0 ipl_height];[transpose(degrees_interp) d_R_interp]]);

% Calculate water vapor partial pressure
e = 6.108 * RH * exp((17.15*T-4684)/(T - 38.45));


% Saastamoinen model
tropoDelaySaas = zeros(size(Zenit_angle));
for i = 1:length(Zenit_angle)
    Zen = Zenit_angle(i);
    tropoDelaySaas(i) = (0.002277/cosd(Zen))*(P+(1255/T+0.05)*e-pow2(tand(Zen)));
end

% corrected model
tropoDelayCorr = zeros(size(Zenit_angle));
for i = 1:length(Zenit_angle)
    Zen = Zenit_angle(i);
    [~, height_idx1] = min(abs(ipl_height - Height));
    [~, height_idx2] = min(abs(B_table.Height - Height));
    [~, zenith_idx] = min(abs(degrees_interp - Zen));
    B_val = B_table.B(height_idx2);
    d_R_val = d_R_interp(zenith_idx,height_idx1);
    tropoDelayCorr(i) = (0.002277/cosd(Zen))*(P+(1255/T+0.05)*e-B_val*pow2(tand(Zen)))+d_R_val;
end

% Return output arguments
delaySaas = tropoDelaySaas;
delayCorr = tropoDelayCorr;



%{
% Plot the interpolated values
plot(height, B, 'o', ipl_height, B_interp);
xlabel('Height');
ylabel('B');
legend('Original', 'Interpolated');
%}

%{
% Plot interpolated data
figure;
pcolor(Xq,Yq,d_R_interp);
shading interp;
ylabel('Degrees');
xlabel('Heights');
colorbar;
%}

%{ 
temp = linspace(-30+273, 30+273, 30); % temperature range from -30 to 30 Celsius
rh = linspace(0, 100, 30); % relative humidity range from 0% to 100%


TD = zeros(length(temp), length(rh)); % initialize tropospheric delay matrix

Zen = 60;
for i = 1:length(temp)
    for j = 1:length(rh)
        e = 6.108 * rh(j)/100 * exp((17.15*temp(i)-4684)/(temp(i) - 38.45));

        [~, height_idx1] = min(abs(ipl_height - Height));
        [~, height_idx2] = min(abs(B_table.Height - Height));
        [~, zenith_idx] = min(abs(degrees_interp - Zen));
        B_val = B_table.B(height_idx2);
        d_R_val = d_R_interp(zenith_idx,height_idx1);
        
        TD(i,j) = (0.002277/cosd(Zen))*(P+(1255/temp(i)+0.05)*e-B_val*pow2(tand(Zen)))+d_R_val
    end
end
surf(temp, rh, TD);
xlabel('Temperature (K)');
ylabel('Relative Humidity (%)');
zlabel('Tropospheric Delay (m)');
title(strcat('Delay for a zenith of  ', string(Zen), ' degrees.'))

colorbar;
xlim([-30+273, 30+273]);
ylim([0, 100]);
%}
end