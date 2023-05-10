clear;clc

c	= 299792458;        %[m/s]
lambda_1  =	0.190293673;    %[m]




% Reciever coordinates
rec_name = {'REF'; 'ROV'};
rec_x  = {3099255449; 3099279.721};
rec_y  = {1014087.319; 1013349.086};
rec_z = {5463148.722; 5463275.806};
rec_clock_corr = {-4.62695e-04; -2.98880e-04};

rec_data = table(rec_name, rec_x , rec_y, rec_z, rec_clock_corr)



% Satellite coordinates
PRN     = {4; 5; 6; 20; 24; 25};
sat_x   = {-9640637.245; 11580980.28; 22299882.17; 15014113.59; 2808519.206; 8650991.184};
sat_y   = {11963653.42; 18202692.29; -922322.5278; 5240279.54; 17209982.45; -12956546.13};
sat_z   = {21764864.54; 15492885.83; 14645540.94; 21124797.43; 20223860.7; 21442358.52};
Ro_dot  = {448.5962979; 479.2977241; -446.2696034; 54.32488024; 50.54129639; -372.5738339};
Ro_REF  = {23408188.07; 21632913.66; 21371250.34; 20127388.69; 21915054.79; 21939371.84};
Ro_ROV  = {23408458.11; 21633431.8; 21371107.05; 20127430.45; 21915515.1; 21938803.05};
P1_ROV  = {23490620.344; 21711251.641; 21416660.203; 20201862.078; 21835232.258; 22030470.547};  
P1_REF  = {23539463.375; 21759846.219;21465916.188; 20250933.188; 21883884.656; 22080151.719};
FI_ROV  = {823151.021; 896009.121; -984361.760; 12789.109; -3418.439; -704176.244};
FI_REF  = {1655955.588; 1790394.946; -2239193.328; -104773.368; -132527.399; -1897452.854};

sat_data = table(PRN, sat_x , sat_y, sat_z, Ro_dot, Ro_REF, Ro_ROV, P1_ROV, P1_REF, FI_ROV, FI_REF) 




P_phase = [ 104166.6667	    -20833.33333	-20833.33333	-20833.33333	-20833.33333;
            -20833.33333	104166.6667	    -20833.33333	-20833.33333	-20833.33333;
            -20833.33333	-20833.33333	104166.6667	    -20833.33333	-20833.33333;
            -20833.33333	-20833.33333	-20833.33333	104166.6667	    -20833.33333;
            -20833.33333	-20833.33333	-20833.33333	-20833.33333	104166.6667]




P_code = [  4.62962963	    -0.925925926	-0.925925926	-0.925925926	-0.925925926;
            -0.925925926	4.62962963	    -0.925925926	-0.925925926	-0.925925926;
            -0.925925926	-0.925925926	4.62962963	    -0.925925926	-0.925925926;
            -0.925925926	-0.925925926	-0.925925926	4.62962963	    -0.925925926;
            -0.925925926	-0.925925926	-0.925925926	-0.925925926	4.62962963]
%parameters 
lamb=lambda_1 %lambda
RoDot=cell2mat(sat_data.Ro_dot)
P1a= cell2mat(sat_data.P1_REF)
ph_sa1= cell2mat(sat_data.FI_REF)
P1b= cell2mat(sat_data.P1_ROV)
ph_sb1= cell2mat(sat_data.FI_ROV)

Xs = cell2mat(sat_data.sat_x)
Ys = cell2mat(sat_data.sat_y)
Zs = cell2mat(sat_data.sat_z)

Xb = cell2mat(rec_data.rec_x(2))
Yb = cell2mat(rec_data.rec_y(2))
Zb = cell2mat(rec_data.rec_z(2))

da=cell2mat(rec_data.rec_clock_corr(1))
db=cell2mat(rec_data.rec_clock_corr(2))

roROV=cell2mat(sat_data.Ro_ROV)
roREF=cell2mat(sat_data.Ro_REF)

P=[[P_phase;zeros(5)] [zeros(5);P_code]]
% lambda 1 from table

% step 1 phase correction and code correction
% rx A code correction
Psa=P1a+RoDot*da
    % read P1 from table REF
    % read RoDot from table
    
    % read Clock corr REF [s] from table
% rx A phase correction
ph_sa=ph_sa1+RoDot*da/lamb
% rx B code correction
Psb=P1b+RoDot*db
    % read P1 from table ROV
    % read RoDot from table (the same)
    
    % read Clock corr ROV [s] from table
% rx A phase correction
ph_sb=ph_sb1+RoDot*db/lamb

% step 2 single difference
Psab=Psb-Psa
ph_sab=ph_sb-ph_sa

%step 3 PRN 20 is our reference satellite and do the double difference
% 4 5 6 24 25
Psab_d=Psab-Psab(4)
ph_sab_d=ph_sab-ph_sab(4)

Psab_d(Psab_d==0) = [];
ph_sab_d(ph_sab_d==0) = [];


%step 4 
ax=-(Xs-Xb)./roROV
ay=-(Ys-Yb)./roROV
az=-(Zs-Zb)./roROV
    %location from satellites and rover
    %roROV [m] from table
ro_ab0=roROV-roREF
%step 5 double difference 
ax_d=ax-ax(4)
ay_d=ay-ay(4)
az_d=az-az(4)

ax_d(ax_d==0)=[];
ay_d(ay_d==0)=[];
az_d(az_d==0)=[];

ro_ab0_d=ro_ab0-ro_ab0(4);

ro_ab0_d(ro_ab0_d==0)=[];


% step 6 build A matrix and vector L


L=[lamb*ph_sab_d-ro_ab0_d ; Psab_d-ro_ab0_d ]


amid=[ax_d ay_d az_d]
A=[[amid;amid] [diag(lamb*[1 1 1 1 1]); zeros(5)]]
X_v=inv(A'*P*A)*A'*P*L

%step 7 Correct the approximate position of the rover
X=Xb+X_v(1)
Y=Yb+X_v(2)
Z=Zb+X_v(3)

%step 8 variance-covariance matrix

Qx=inv(A'*P*A) %Qx is 8*8 matrix, to get ambiguity fetch 5*5 from the rest part of Qx
Qahat=Qx([4:8],[4:8])
ahat = X_v(end-4:end)

%step 9 use the LAMBDA function
[afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu]=LAMBDA(ahat,Qahat);

%step 10 




