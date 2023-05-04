clear
clc
% read satellite
M=readtable("output.csv");
D=readtable("delays.csv");

filename = 'rinex.xlsx';
table_data = readtable(filename);
% step 9: troposphere correction
T=D.d_trop_m_; %m
% step 10: ionosphere correction
I=D.d_ion_m_; %m
% step 11 and 16: 
e=1e-5; %accept error level
err_pre=0;
con=1;
cnt=0;
dta=0;
%Parameters 
Omega_dot=7.2921151467e-5; %[m^3/s^2] earth rotation rate
c=299792458; %[m/s] speed of light
%X, Y, Z for receiver status
X=3104219.4530;
Y=998383.9820;
Z=5463290.5080;
Xa=X;
Ya=Y;
Za=Z;
for i = 1:height(table_data)-1
    Psa(i)=table_data.Pa(i);
end
Psa=Psa';

while(con>e) 
    d_tsa=Psa/c;  %time difference between S and R
    xa=Xa-Omega_dot*Ya*d_tsa;
    ya=Ya+Omega_dot*Xa*d_tsa;
    za=Za;
    %satellite position
    Xs=M.X;
    Ys=M.Y;
    Zs=M.Z;
    d_tL1=M.dts_Li;
    
    mx=Xs-xa;
    my=Ys-ya;
    mz=Zs-za;
    
    dist_sa=sqrt(mx.^2+my.^2+mz.^2) % convert into matrix calculation
    %first time neglect?
    
    
    %%%%%%%%%%%%%%%%% <= todoooooo

    % step 12 Compute L vector eq 19
  
    L=Psa-dist_sa+c*d_tL1-I-T;
	    %d_tL1: satellite clock error
	    %I: Ionosphere delay
	    %Troposphere delay
  
    % step 13 Compute A matrix eq20 eq12
    A=ones(12,4);
	%eq 12
    ax=-(Xs-Xa)./dist_sa;
    ay=-(Ys-Xa)./dist_sa;
    az=-(Zs-Xa)./dist_sa;
	%eq20
    A(:,1)=ax;
    A(:,2)=ay;
    A(:,3)=az;
    % step 14 unknown parameters in eq18
    Qx=inv(A'*A);
    X_v=Qx*A'*L;
    v=A*X_v-L;
    err=v'*v;
    con=abs(err-err_pre);
    err_pre=err;
    %some quallity things
    s0=sqrt(err/(12-4));
    sx=s0*sqrt(Qx(1,1));
    sy=s0*sqrt(Qx(2,2));
    sz=s0*sqrt(Qx(3,3));
    sc=s0*sqrt(Qx(4,4));
    PDOP=sqrt(Qx(1,1)+Qx(2,2)+Qx(3,3));
    % step 15 update the receiver coord. eq22
    Xa=Xa+X_v(1);
    Ya=Ya+X_v(2);
    Za=Za+X_v(3);
    dta=X_v(4)/c;
    % step 16 iteration
    %Psa=dist_sa+c*dta-c*d_tL1+I+T; %correct till here
    cnt=cnt+1;
end

Xa
Ya
Za
sx
sy
sz
sc
dta
PDOP
