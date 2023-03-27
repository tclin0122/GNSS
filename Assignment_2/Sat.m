%% Parameters
% Global
c = 299792458;               %[m/s]
pi = 3.1415926535898;
mu = 3.986005e14;            %[m3/s2]
Omega_E = 7.2921151467e-5;   %[rad/s]

% Step 3
a_f0 = 1.093372702600e-06; %1-2
a_f1 = 2.955857780760e-12;  %1-3
a_f2 = 000000000000e+00;    %1-4
T_GD = -9.313225746150e-10; %7-3

% Step 5

A_sqrt = 5.153573421480e+03; %3-4
A =  A_sqrt^2;
t_oe = 9.360000000000e+04;  %4-1

% Step 6

e = 9.439246961850e-03; %3-2
F = -4.442807633e-10; %[s /m^0.5]
dn = 4.158030341510e-09; %2-3
M_0 = -1.234439725290e+00; %2-4



%step 7
%line 10
w=-1.409026029280e+00; %5-3
%line 11
Cus=7.957220077510e-06; %3-3
Cuc=5.569308996200e-06; %3-1
%line 12
Crs=1.055312500000e+02; %2-2
Crc=2.370625000000e+02; %5-2
%line 13
Cis=1.303851604460e-07; %4-4
Cic=-1.676380634310e-07; %4-2
%line 16
i0=9.719990595480e-01; %5-1
IDOT=-1.821504444400e-11;%6-1
%line 18
Omega_0=2.953800868980e+00;  %4-3
Omega_dot=-7.928544541420e-09; %5-4
Omega_dot_e=7.2921151467e-5; % constant




%% 1. Compute signal propagation time

tt_A = 1 * 24 * 3600 + 3600;                   % passed time since GPS week started [s]
Ps_Att = 25001256.67943;                       % [km]

dts_A = Ps_Att/c


%% 2. Compute signal transmission time

tts = tt_A - Ps_Att / c

%% 3. Compute satellite clock correction s L1 t by (27) and (28), neglect Δtr.

dt_r = 0;                                       % neglectd for now

t_oc = 86400 + 3600 * 2;
dt_sv = a_f0 + a_f1*(tts-t_oc) + a_f2*(tts-t_oc)^2 + dt_r
dts_L1 = dt_sv - T_GD


%% 4. Compute ts by (15) using the correction from the step 3.

ts = tts - dts_L1

%% 5. Compute eccentric anomaly (Table 2 - line 4-8)

% Computed mean motion
n_0 = sqrt(mu/A^3)

t_k = ts - t_oe

%% 6. Compute Δtr by (29), SVt by (28), s L1 t by (27) and ts by (15).


n = n_0 + dn

M_k = M_0 + n * t_k
E_0 = M_k

E_ip = E_0
for i = 10                      % TODO proper iteration number
   E_i = M_k + e * sin(E_ip) 
   E_ip = E_i
end

E_k = E_i;
dtr = F * e * A_sqrt * sin(E_k)

% for for dtr

%% 7. Compute satellite coordinates Xs, Ys, Zs, for time s t - Table 2 (line 4-19). Update the eccentric anomaly computed in step 5.

%% Equations for coordinates output
%line 9
a=(sqrt(1-e^2)*sin(Ek))/(1-e*cos(Ek));
b=(cos(Ek)-e)/(1-e*cos(Ek));
vk=atan2(a,b) %recommand to use atan

%line 10
latk=vk+w

%line 11
d_uk=Cus*sin(2*latk)+Cuc*cos(2*latk)
%line 12
d_rk=Crs*sin(2*latk)+Crc*cos(2*latk)
%line 13
d_ik=Cis*sin(2*latk)+Cic*cos(2*latk)

%line 14
uk=latk+d_uk
%line 15
rk=A*(1-e*cos(Ek))+d_rk
%line 16
ik=i0+d_ik+IDOT*tk
%line 17
xkprime=rk*cos(uk)
ykprime=rk*sin(uk)
%line 18
Omega_k=Omega_0+(Omega_dot-Omega_dot_e)*tk-Omega_dot_e*toe
%line 19
xk=xkprime*cos(Omega_k)-ykprime*cos(ik)*sin(Omega_k)
yk=xkprime*sin(Omega_k)+ykprime*cos(ik)*cos(Omega_k)
zk=ykprime*sin(ik)