clear clc


%% Positon data
filename_pos = 'tester.csv';
fileID = fopen(filename_pos,'r');
headerLine = fgetl(fileID);
columnNames = strsplit(headerLine, ',');
Pos_data = readtable(filename_pos);
fclose(fileID);

%% data from website
TECU=5.3;


d_trop_func = getTroposfericDelay(60,3.2)

d_iono_func = getIonosphericDelay(TECU,61.3512)