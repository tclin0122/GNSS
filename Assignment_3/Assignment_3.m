clear clc


%% Positon data
filename_pos = 'output.csv';
fileID = fopen(filename_pos,'r');
headerLine = fgetl(fileID);
columnNames = strsplit(headerLine, ',');
Pos_data = readtable(filename_pos);
fclose(fileID);
%% Readfile 

%Receiver position
x=3104219.4530;
y=998383.9820;
z=5463290.5080;
[lat, lon, h]=Ellop2Car(x, y, z)
Xr=[x,y,z];
[Zenit_angle]=relative_pos(x,y,z,lat,lon)


%% data from website
TECU=5.3;


d_trop_func = getTroposfericDelay(61.25,3.2)

d_iono_func = getIonosphericDelay(TECU,Zenit_angle)

