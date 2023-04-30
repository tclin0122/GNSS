clear clc


%% Positon data
filename_pos = 'tester.csv';
fileID = fopen(filename_pos,'r');
headerLine = fgetl(fileID);
columnNames = strsplit(headerLine, ',');
Pos_data = readtable(filename_pos);
fclose(fileID);




d_trop_func = getTroposfericDelay(60,3.2)
