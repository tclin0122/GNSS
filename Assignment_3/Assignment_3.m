clear clc

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
B_interp = interp1(height, B, new_height); % interpolate B values at new height values

% Plot the interpolated values
%plot(height, B, 'o', new_height, B_interp);
%xlabel('Height');
%ylabel('B');
%legend('Original', 'Interpolated');

B_table = array2table([transpose(new_height) transpose(B_interp)]);
B_table.Properties.VariableNames(1:2) = {'Height','B'};


%% d_R data

d_R_data = readmatrix('d_R_data.csv'); % Read CSV file
heights = d_R_data(1,2:end); % Extract heights
degrees = d_R_data(2:end,1); % Extract degrees
values = d_R_data(2:end,2:end); % Extract data values

[X,Y] = meshgrid(heights,degrees); % Create grid of original data points
degrees_interp = min(degrees):0.5:max(degrees);
heights_interp = min(heights):0.5:max(heights);
[Xq,Yq] = meshgrid(heights_interp,degrees_interp); % Create query grid

d_R_interp = interp2(X,Y,values,Xq,Yq);

% Plot interpolated data
%figure;
%pcolor(Xq,Yq,d_R_interp);
%shading interp;
%ylabel('Degrees');
%xlabel('Heights');



d_R_table = array2table([[0 heights_interp];[transpose(degrees_interp) d_R_interp]]);
%d_R_table.Properties.VariableNames(1:2) = {'Height','B'}


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
Hoff_height = 4.99999;  % TODOOOOO 
Hoff_zenith = 60;
B_val = B_table.B(B_table.Height == Hoff_height)
d_R_val = d_R_interp(find(degrees_interp > Hoff_zenith && degrees_interp <= Hoff_zenith),find(heights_interp == Hoff_height))
%d_R_val = 0*d_R_interp;
d_trop = (0.002277/cosd(Zen))*(P+(1255/T+0.05)*e-B_val*pow2(tand(Zen)))+d_R_val



%elev to Zenit 

