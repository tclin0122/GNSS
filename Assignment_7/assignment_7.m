clc 
clear all
%% Assignment 7 Kalman filter
M=readtable("Measeurements.xlsx");
l=length(M.Var1); % will be our loop interation time
% Kalman filter matrix
F=zeros(4);
F(1,3)=1;
F(2,4)=1;

G=zeros(4,2);
G(3,2)=1;
G(4,1)=1;

%Standard deviation of measured coordinates 3 m/s
% Standard deviation of measured abs. velocity 0.5 m/s
r1=3^2
r3=0.5^2
R=diag([r1 r1 r3])

dt=2 %second
I=diag([1 1 1 1]);
T=I+F*dt;

%initial value
% Standard deviation of initial velocity 3 m/s
% Standard deviation of initial coordinates 10 m
x0=[M.Var2(1); M.Var3(1); 3.53; 0.86];
Q0=I;
Q0(1,1)=10^2;
Q0(2,2)=10^2;
Q0(3,3)=3^2;
Q0(4,4)=3^2;

%Covariance matrix of model%s noise Qk
%we use PSD = 0.01 m^2 s^-3
qe=0.01;
qn=0.01;
Qk=[qe*dt^3/3 0 qe*dt^2/2 0;0 qn*dt^3/3 0 qn*dt^2/2; qe*dt^2/2 0 qe*dt 0; 0 qn*dt^2/2 0 qn*dt];


% step 1 Initialization
x0
Q0
x_f=zeros(l,4); % store data
x_f(1,:)=x0'
Q_f=zeros(4*25,4);% store data Q0*25
Q_f(1:4,1:4)=Q0
Q_p=zeros(4*25,4);% store data Q0*25
Q_p(1:4,1:4)=Q0
for i=1:(l-1)
% step 2 Time propagation <= start the loop from here
x_p=T*x0
v1=sqrt(x_p(3)^2+x_p(4)^2)
Qx1=T*Q0*T'+Qk
Q_p(i*4+1:i*4+4,1:4)=Qx1

H=[1 0 0 0; 0 1 0 0; 0 0 x_p(3)/v1 x_p(4)/v1]
% step 3 Gain calculation
K1=Qx1*H'*inv(R+H*Qx1*H')

% step 4 Measurement update
L1=[M.Var2(i+1); M.Var3(i+1); M.Var4(i+1)]
hkxk_p=[x_p(1);x_p(2);v1]
x_p=x_p+K1*[L1-hkxk_p]
x0=x_p
x_f(i+1,:)=x_p'
% Step 5 Covariance update
Q0=[I-K1*H]*Qx1

Q_f(i*4+1:i*4+4,1:4)=Q0
end

%plotting the result

% success without looping
% read file
tr=readtable("True_value.xlsx")
figure(1)
plot(x_f(:,1),x_f(:,2))
hold on
plot(tr.e_m_,tr.n_m_)
legend('filter','observation')
hold off

% difference
for i=1:(l)
    e_d(i)=x_f(i,1)-tr.e_m_(i);
    n_d(i)=x_f(i,2)-tr.n_m_(i);
    ve_d(i)=x_f(i,3)-tr.ve_m_s_(i);
    vn_d(i)=x_f(i,4)-tr.vn_m_s_(i);
end
figure(2)
plot(tr.Time_s_,e_d')
hold on
plot(tr.Time_s_,n_d')
plot(tr.Time_s_,ve_d')
plot(tr.Time_s_,vn_d')
legend('e_d (m)','n_d (m)','ve_d (m/s)','vn_d (m/s)')
hold off


%% Smoothing <=TOOOOOOODOOOOOOO
x_s=zeros(l,4);
x_s(l,:)=x_f(l,:);
a=l;
x_f(a,:)
while(a>=2)
    D=Q_f((a-1)*4-3:(a-1)*4,1:4)*T'*inv(Q_p((a)*4-3:(a)*4,1:4));
    Q_f((a-1)*4-3:(a-1)*4,1:4)
    x_s(a-1,:)=(x_f(a-1,:)'+D*((x_s(a,:)-x_f(a,:))'))';
    %x_s(i-1,:)=
    a=a-1;
end

%plotting the result

% success without looping
% read file

figure(3)
plot(x_s(:,1),x_s(:,2))
hold on

plot(tr.e_m_,tr.n_m_)
legend('smooth','filter','observation')
hold off

% difference
T=readtable("True_value.xlsx")
for i=1:(l)
    se_d(i)=x_s(i,1)-tr.e_m_(i);
    sn_d(i)=x_s(i,2)-tr.n_m_(i);
    sve_d(i)=x_s(i,3)-tr.ve_m_s_(i);
    svn_d(i)=x_s(i,4)-tr.vn_m_s_(i);
end

figure(4)
plot(tr.Time_s_,se_d')
hold on
plot(tr.Time_s_,sn_d')
plot(tr.Time_s_,sve_d')
plot(tr.Time_s_,svn_d')
legend('e_d (m)','n_d (m)','ve_d (m/s)','vn_d (m/s)')
hold off

