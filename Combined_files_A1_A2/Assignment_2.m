clear all

%% Parameters
% Global
c = 299792458;               %[m/s]
pi = 3.1415926535898;
mu = 3.986005e14;            %[m3/s2]
Omega_E = 7.2921151467e-5;   %[rad/s]


%% File reading
filename = 'rinex.xlsx';
table_data = readtable(filename);
disp(table_data);

%% File writting
output_filename = 'output.csv';
output_file = fopen(output_filename, 'w');
fprintf(output_file, 'PRN,X,Y,Z,dts_Li\n');


% Iterate through each row of the table
for i = 1:height(table_data)-1
    % Get the data from the i-th row of the table
    row_data = table_data(i,:);
    
    % Process the row_data however you want
    % For example, you can display the values in the row
    %disp(row_data);


    satnum = table_data.PRN(i);             % Satellite number
    
    % Step 3
    a_f0 = table_data.a_f0(i);              %1-2
    a_f1 = table_data.a_f1(i);              %1-3
    a_f2 = table_data.a_f2(i);              %1-4
    T_GD = table_data.T_GD(i);              %7-3
    
    % Step 5    
    A_sqrt = table_data.sqrt_A(i);          %3-4
    A =  A_sqrt^2;
    t_oe = table_data.t_oe(i);              %4-1
    
    % Step 6    
    e = table_data.e(i);                    %3-2
    F = -4.442807633e-10;                   %[s /m^0.5]
    dn = table_data.dn(i);                  %2-3
    M_0 = table_data.M_0(i);                %2-4
    
    
    
    %step 7
    %line 10
    w=table_data.w(i); %5-3
    %line 11
    Cus=table_data.Cus(i); %3-3
    Cuc=table_data.Cuc(i); %3-1
    %line 12
    Crs=table_data.Crs(i); %2-2
    Crc=table_data.Crc(i); %5-2
    %line 13
    Cis=table_data.Cis(i); %4-4
    Cic=table_data.Cic(i); %4-2
    %line 16
    i0=table_data.i0(i); %5-1
    IDOT=table_data.IDOT(i);%6-1
    %line 18
    Omega_0=table_data.Omega_0(i);  %4-3
    Omega_dot=table_data.Omega_dot(i); %5-4
    Omega_dot_e=7.2921151467e-5; % constant
    
    
    
    
    %% 1. Compute signal propagation time
    
    tt_A = 1 * 24 * 3600 + 3600;                   % passed time since GPS week started [s]
    Ps_Att = 25001256.67943;                       % [km]
    
    dts_A = Ps_Att/c;
    
    
    %% 2. Compute signal transmission time
    
    tts = tt_A - Ps_Att / c;
    
    %% 3. Compute satellite clock correction s L1 t by (27) and (28), neglect Δtr.
    
    dt_r = 0;                                       % neglectd for now
    
    
    for a = 1:20 
        t_oc = 86400 + 3600 * 2;
        dt_sv = a_f0 + a_f1*(tts-t_oc) + a_f2*(tts-t_oc)^2 + dt_r;
        dts_L1 = dt_sv - T_GD;
        
        
        %% 4. Compute ts by (15) using the correction from the step 3.
        
        ts = tts - dts_L1;
        
        %% 5. Compute eccentric anomaly (Table 2 - line 4-8)
        
        % Computed mean motion
        n_0 = sqrt(mu/A^3);
        t_k = ts - t_oe;
        
        %% 6. Compute Δtr by (29), SVt by (28), s L1 t by (27) and ts by (15).
        
        n = n_0 + dn;
        
        M_k = M_0 + n * t_k;
        E_0 = M_k;
        
        E_ip = 0;
        E_i = E_0;
        %cp=abs(E_i-E_ip)
        
        for a = 1:20 
	        E_i=M_k+e*sin(E_ip);
	        %cp=abs(E_i-E_ip)
            E_ip = E_i;
        end
        
        
        
        E_k = E_i;
        dt_r = F * e * A_sqrt * sin(E_k);
    end
    
    
    %% 7. Compute satellite coordinates Xs, Ys, Zs, for time s t - Table 2 (line 4-19). Update the eccentric anomaly computed in step 5.
    
    %% Equations for coordinates output
    %line 9
    a=(sqrt(1-e^2)*sin(E_k))/(1-e*cos(E_k));
    b=(cos(E_k)-e)/(1-e*cos(E_k));
    vk=atan2(a,b); %recommand to use atan
    
    %line 10
    latk=vk+w;
    
    %line 11
    d_uk=Cus*sin(2*latk)+Cuc*cos(2*latk);
    %line 12
    d_rk=Crs*sin(2*latk)+Crc*cos(2*latk);
    %line 13
    d_ik=Cis*sin(2*latk)+Cic*cos(2*latk);
    
    %line 14
    uk=latk+d_uk;
    %line 15
    rk=A*(1-e*cos(E_k))+d_rk;
    %line 16
    ik=i0+d_ik+IDOT*t_k; 
    %line 17
    xkprime=rk*cos(uk);
    ykprime=rk*sin(uk);
    %line 18
    Omega_k=Omega_0+(Omega_dot-Omega_dot_e)*t_k-Omega_dot_e*t_oe; 
    %line 19
    xk=xkprime*cos(Omega_k)-ykprime*cos(ik)*sin(Omega_k);
    yk=xkprime*sin(Omega_k)+ykprime*cos(ik)*cos(Omega_k);
    zk=ykprime*sin(ik);
    
    
    
    fprintf('Result for satellite number ')
    fprintf('%d',satnum)
    fprintf(':')
    fprintf('\n')
    
    fprintf('X= ')
    fprintf('%d', xk)
    fprintf('  Y= ')
    fprintf('%d', yk)
    fprintf('  Z= ')
    fprintf('%d', zk)
    fprintf('\n')
    
    fprintf("Satellite clock correction: ")
    fprintf('%d', dts_L1)
    fprintf('\n')

    fprintf(output_file, '%d,%f,%f,%f,%f\n', satnum, xk, yk, zk, dts_L1)

end

% Close the output text file
fclose(output_file);