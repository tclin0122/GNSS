function [outputArg1,outputArg2] = getTroposfericDelay(Zenit_angle,Height)

%% Parameters
P = 1013;           % [mBar,hPa] 
T = 273.15+ 18;     % [Kelvin]
RH = 0.5;

P = 1013.25;           % [mBar,hPa] 
T = 273.15;     % [Kelvin]
RH = 0;
Zen = 0;

%% Simulation parameters
interp_resolution = 0.05;    

%% File handling
B_data = readtable('B_data.csv');
d_R_data = readmatrix('d_R_data.csv');






%% B data
% Extract the height and B columns
height = B_data.Height;
B = B_data.B;

old_height = d_R_data(1,2:end); % Extract heights
ipl_height = min(old_height):interp_resolution:max(old_height);
B_interp = interp1(height, B, ipl_height); % interpolate B values at new height values

% Plot the interpolated values
%plot(height, B, 'o', ipl_height, B_interp);
%xlabel('Height');
%ylabel('B');
%legend('Original', 'Interpolated');

B_table = array2table([transpose(ipl_height) transpose(B_interp)]);
B_table.Properties.VariableNames(1:2) = {'Height','B'};


%% d_R data
degrees = d_R_data(2:end,1); % Extract degrees
values = d_R_data(2:end,2:end); % Extract data values

[X,Y] = meshgrid(old_height,degrees); % Create grid of original data points
degrees_interp = min(degrees):0.5:max(degrees)
[Xq,Yq] = meshgrid(ipl_height,degrees_interp); % Create query grid
d_R_interp = interp2(X,Y,values,Xq,Yq);

% Plot interpolated data
%figure;
%pcolor(Xq,Yq,d_R_interp);
%shading interp;
%ylabel('Degrees');
%xlabel('Heights');



d_R_table = array2table([[0 ipl_height];[transpose(degrees_interp) d_R_interp]]);




e = 6.108 * RH * exp((17.15*T-4684)/(T - 38.45));

%Saastamoinen model
troposfericDelaySaas = (0.002277/cosd(Zen))*(P+(1255/T+0.05)*e-pow2(tand(Zen)));


%Hoffman model

% Find the index of the element(s) in Heights with the minimum absolute difference
[~, height_idx1] = min(abs(ipl_height - Height));
[~, height_idx2] = min(abs(B_table.Height - Height));
[~, zenith_idx] = min(abs(degrees_interp - Zenit_angle));



B_val = B_table.B(height_idx2);
d_R_val = d_R_interp(zenith_idx,height_idx1);

troposfericDelayHoff = (0.002277/cosd(Zen))*(P+(1255/T+0.05)*e-B_val*pow2(tand(Zen)))+d_R_val



outputArg1 = troposfericDelaySaas;
outputArg2 = troposfericDelayHoff;


end