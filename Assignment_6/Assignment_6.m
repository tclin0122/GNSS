load("data.m")
%parameters 
lamb=lambda_1 %lambda
RoDot=sat_data.Ro_dot
P1a= sat_data.P1_REF
ph_sa1= sat_data.FI_REF
P1b= sat_data.P1_ROV
ph_sb1= sat_data.FI_ROV

da=rec_data.rec_clock_corr(1)
db=rec_data.rec_clock_corr(2)

roROV=sat_data.Ro_ROV
roREF=sat_data.Ro_REF

P=[[P_phase;zeros(5)] [zeros(5);P_code]]
% lambda 1 from table

% step 1 phase correction and code correction
% rx A code correction
Psa=P1a+RoDot*da
    % read P1 from table REF
    % read RoDot from table
    
    % read Clock corr REF [s] from table
% rx A phase correction
ph_sa=ph_sa1+RoDot*da/lam
% rx B code correction
Psb=P1b+RoDot*db
    % read P1 from table ROV
    % read RoDot from table (the same)
    
    % read Clock corr ROV [s] from table
% rx A phase correction
ph_sb=ph_sb1+RoDot*db/lam

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
ax=-(Xs-Xb)/roROV
ay=-(Ys-Yb)/roROV
az=-(Zs-Zb)/roROV
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


L=[lam*ph_sab_d-ro_ab0_d ; Psab_d-ro_ab0_d ]


amid=[ax ay az]
A=[[amid;amid] [diag(lam*[1 1 1 1 1]); zeros(5)]]
X_v=inv(A'*P*A)*A'*P*L

%step 7 Correct the approximate position of the rover
X=Xb+X_v(1)
Y=Yb+X_v(2)
Z=Zb+X_v(3)

%step 8 variance-covariance matrix

Qahat=inv(A'*P*A) %Qx is 8*8 matrix, to get ambiguity fetch 5*5 from the rest part of Qx
ahat=Qx([4:8],[4:8])

%step 9 use the LAMBDA function
[afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu]=LAMBDA(ahat,Qahat);

%step 10 




