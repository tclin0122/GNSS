function [zij] = relative_pos(x,y,z,lat,lon)
%RELATIVE_POS Summary of this function goes here
%   Detailed explanation goes here
M=readtable("output.csv");
Xij=[M.X-x,M.Y-y,M.Z-z];

R=[-sind(lat)*cosd(lon), -sind(lon), cosd(lat)*cosd(lon);-sind(lat)*sind(lon),cosd(lon),cosd(lat)*sind(lon);cosd(lat),0,sind(lat)];
sij=zeros(1,29);
alpha_ij=zeros(1,12);
zij=zeros(1,12);
xij=zeros(3,12);
for i=1:12
    xij(:,i)=R'*(Xij(i,:)');
    sij(i)=(sqrt(power(xij(1,i), 2) + power(xij(2,i), 2) + power(xij(3,i), 2)));  % spatial distance
    alpha_ij(i)=(rad2deg(atan2(xij(2,i), xij(1,i))));  % azimuth
    zij(i)=(rad2deg(acos(xij(3,i) / (sqrt(power(xij(1,i), 2) + power(xij(2,i), 2) + power(xij(3,i), 2))))));
end
end

