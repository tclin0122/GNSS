
%parameters 
lam %lambda
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

%step 5 double difference 
ax_d=ax-ax(4)
ay_d=ay-ay(4)
az_d=az-az(4)





