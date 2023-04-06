% step 9: troposphere correction
T=2.39 %m
% step 10: ionosphere correction
I=2 %m
% step 11: 
%Parameters 
Omega_dot=7.2921151467e-5 %[m^3/s^2] earth rotation rate
c=299792458 %[m/s] speed of light
d_tsa=(ta-ts);  %time difference between S and R
%X, Y, Z for receiver status
X=3104219.4530
Y=998383.9820
Z=5463290.5080

xa=X-Omega_dot*Y*d_tsa
ya=Y-Omega_dot*X*d_tsa
za=Z

dist_sa=sqrt(pow2(Xs-xa)+pow2(Ys-ya)+pow2(Zs-za)) % convert into matrix calculation

% step 12 Compute L vector eq 19
  L=zeros(1,12)
  for i=(1:12)
  	L(i)=dist_sa(i)-c*d_tsa(i)+c*d_tL1(i)-I-T
  	%d_tL1: satellite clock error
  	%I: Ionosphere delay
  	%Troposphere delay
  end
% step 13 Compute A matrix eq20 eq12
A=ones(12,4)
for i=
	%eq 12
	ax(i)=-(Xs(i)-X)/dist_sa(i)
	ay(i)=-(Ys(i)-Y)/dist_sa(i)
	ay(i)=-(Zs(i)-Z)/dist_sa(i)
end
	%eq20
	A(:,1)=ax
	A(:,2)=ay
	A(:,3)=az
% step 14 unknown parameters in eq18
Qx=inv(A'*A)
X=Q*A'*L

% step 15 update the receiver coord. eq22
xa=xa+X(1)
ya=ya+X(2)
za=za+X(3)
% step 16 iteration
