clear all

%% Positon data
filename_pos = 'tester.csv';
fileID = fopen(filename_pos,'r');
headerLine = fgetl(fileID);
columnNames = strsplit(headerLine, ',');
Pos_data = readtable(filename_pos);
fclose(fileID);


%% B data

% Read the CSV file into a table
B_data = readtable('B_data.csv');

% Extract the height and B columns
height = B_data.Height;
B = B_data.B;

% Interpolate the values
new_height = linspace(min(height), max(height), 100); % new height values to interpolate
new_B = interp1(height, B, new_height); % interpolate B values at new height values

% Plot the interpolated values
plot(height, B, 'o', new_height, new_B);
xlabel('Height');
ylabel('B');
legend('Original', 'Interpolated');


%% d_R data

data = readmatrix('d_R_data.csv'); % Read CSV file
degrees = data(1,2:end); % Extract degrees
heights = data(2:end,1); % Extract heights
values = data(2:end,2:end); % Extract data values

[X,Y] = meshgrid(degrees,heights); % Create grid of original data points
[Xq,Yq] = meshgrid(min(degrees):0.5:max(degrees),min(heights):0.5:max(heights)); % Create query grid

interp_values = interp2(X,Y,values,Xq,Yq);

% Plot interpolated data
figure;
pcolor(Xq,Yq,interp_values);
shading interp;
xlabel('Degrees');
ylabel('Heights');




%% Parameters
P = 1013;           % [mBar,hPa] 
T = 273.15+ 18;     % [Kelvin]
RH = 0.5;

P = 1013.25;           % [mBar,hPa] 
T = 273.15;     % [Kelvin]
RH = 0;
Zen = 0;

e = 6.108 * RH * exp((17.15*T-4684)/(T - 38.45));

%Saastamoinen model
d_trop = (0.002277/cosd(Zen))*(P+(1255/T+0.05)*e-pow2(tand(Zen)))

%Hoffman model
%d_trop = (0.002277/cosd(Zen))*(P+(1255/T+0.05)*e-B*pow2(tand(Zen)))+d_R



%elev to Zenit 

