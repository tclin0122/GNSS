%% function
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

