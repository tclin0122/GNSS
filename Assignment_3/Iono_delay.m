clc
clear

M=readtable("output.csv")


%parameter
f=1575.42*1e6;
Re=6371000; %m
hi=350000; %m

%Receiver position
x=3104219.4530
y=998383.9820
z=5463290.5080
Xr=[x,y,z]
Xij=[M.X-x,M.Y-y,M.Z-z]

lat=59.3497
lon=18.0694
R=[-sind(lat)*cosd(lon), -sind(lon), cosd(lat)*cosd(lon);-sind(lat)*sind(lon),cosd(lon),cosd(lat)*sind(lon);cosd(lat),0,sind(lat)]
sij=zeros(1,29);
alpha_ij=zeros(1,12);
zij=zeros(1,12);
xij=zeros(3,12)
for i=1:12
    xij(:,i)=R'*(Xij(i,:)');
    sij(i)=(sqrt(power(xij(1,i), 2) + power(xij(2,i), 2) + power(xij(3,i), 2)));  % spatial distance
    alpha_ij(i)=(rad2deg(atan2(xij(2,i), xij(1,i))));  % azimuth
    zij(i)=(rad2deg(acos(xij(3,i) / (sqrt(power(xij(1,i), 2) + power(xij(2,i), 2) + power(xij(3,i), 2))))));
end

el=90
%1 TEC Unit TECU = 10^16 electrons/m²
TECU=5.3;
TEC=TECU*10^16;
zen=90-el;
OF=(1-((Re*sind(zij)/(Re+hi)).^2)).^(-0.5)
dg=(40.3/(f^2))*TEC*OF

Res=zeros(3,12)
Res(1,:)=M.PRN
Res(2,:)=zij
Res(3,:)=dg




function [lat0, lon, h]=Ellop2Car(X, Y, Z)
    a = 6378137;  % semi major axis
    fe = 1 / 298.257222101;
    b = a - fe * a; % semi minor axis
    e2 = (power(a, 2) - power(b, 2)) / power(a, 2);
    p = sqrt(power(X, 2) + power(Y, 2));
    lat0 = rad2deg(atan2((Z / p) * (power((1 - e2), -1)),1));
    lat = 0;
    while (lat~=lat0)
        N0 = power(a, 2) / sqrt(power(a, 2) * power(cosd(lat0), 2) + power(b, 2) * power(sind(lat0), 2));
        h = p / cosd(lat0) - N0;
        lat = rad2deg(atan2((Z / p) * power(1 - e2 * (N0 / (N0 + h)), -1),1));
        lat0 = lat;
        lon = rad2deg(atan2(Y,X));
        h = p / cosd(lat0) - N0;
    end

end