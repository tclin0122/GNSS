import numpy as np
import pandas as pd
import array as arr

##parameters
# latitude - phi
LaD = 59  # 47
LaM = 20  # 0
LaS = 59  # 0
# longitude - lambda
LoD = 18  # 15
LoM = 4  # 0
LoS = 10  # 0

lat = LaD + LaM / 60 + LaS / 3600
lon = LoD + LoM / 60 + LoS / 3600
h = 62 + 23 + 100
print(lat, lon, h)


def Car2Ellip(lat, lon, h):
    a = 6378137  # semi major axis
    fe = 1 / 298.257222101
    b = a - fe * a  # semi minor axis
    N = pow(a, 2) / np.sqrt(pow(a, 2) * pow(np.cos(np.deg2rad(lat)), 2) + pow(b, 2) * pow(np.sin(np.deg2rad(lat)), 2))
    X = (N + h) * np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(lon))
    Y = (N + h) * np.cos(np.deg2rad(lat)) * np.sin(np.deg2rad(lon))
    Z = (pow(b / a, 2) * N + h) * np.sin(np.deg2rad(lat))
    # print(X,Y,Z)
    return X, Y, Z


[X, Y, Z] = Car2Ellip(lat, lon, h)

# Part B - From global to local
# use the xlsx file
# import the GPS satellite position
# need to use the right place
dataframe1 = pd.read_excel('./GPS_position.xlsx', index_col=None, header=None)
Xs = pd.read_excel('./GPS_position.xlsx', usecols='B', index_col=None, header=None).values
Ys = pd.read_excel('./GPS_position.xlsx', usecols='C', index_col=None, header=None).values
Zs = pd.read_excel('./GPS_position.xlsx', usecols='D', index_col=None, header=None).values
print(dataframe1)
print(Xs[0])
# define another position to test the code 100m above the point

[X, Y, Z] = Car2Ellip(lat, lon, h)


# print(X,Y,Z)
# invert cart. back to ellipse

def Ellop2Car(X, Y, Z):
    a = 6378137  # semi major axis
    fe = 1 / 298.257222101
    b = a - fe * a  # semi minor axis
    e2 = (pow(a, 2) - pow(b, 2)) / pow(a, 2)
    p = np.sqrt(pow(X, 2) + pow(Y, 2))
    lat0 = np.rad2deg(np.arctan2((Z / p) * (pow((1 - e2), -1)), 1))
    lat = 0
    while (lat != lat0):
        N0 = pow(a, 2) / np.sqrt(
            pow(a, 2) * pow(np.cos(np.deg2rad(lat0)), 2) + pow(b, 2) * pow(np.sin(np.deg2rad(lat0)), 2))
        h = p / np.cos(np.deg2rad(lat0)) - N0
        lat = np.rad2deg(np.arctan2((Z / p) * pow(1 - e2 * (N0 / (N0 + h)), -1), 1))
        lat0 = lat
    lon = np.rad2deg(np.arctan2(Y, X))
    h = p / np.cos(np.deg2rad(lat0)) - N0
    return lat0, lon, h


[lat, lon, h] = Ellop2Car(X, Y, Z)
print(lat, lon, h)
# start from global to local part
Xi = np.array([X, Y, Z])

Xj = np.array([Xs * 1000, Ys * 1000, Zs * 1000])

ni = np.array([-np.sin(np.deg2rad(lat)) * np.cos(np.deg2rad(lon)), -np.sin(np.deg2rad(lat)) * np.sin(np.deg2rad(lon)),
               np.cos(np.deg2rad(lat))])
ei = np.array([-np.sin(np.deg2rad(lon)), np.cos(np.deg2rad(lon)), 0])
ui = np.array([np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(lon)), np.cos(np.deg2rad(lat)) * np.sin(np.deg2rad(lon)),
               np.sin(np.deg2rad(lat))])
R = np.array([ni, ei, ui])
print(R)
cnt = 0
for i in range(0, 29):
    Xij = np.array([Xj[0][i] - X, Xj[1][i] - Y, Xj[2][i] - Z])
    xij = np.dot(np.transpose(R), Xij)
    # initialize dynamic array
    sij = []
    alpha_ij = []
    zij = []
    sij.append(np.sqrt(pow(xij[0], 2) + pow(xij[1], 2) + pow(xij[2], 2)))  # spatial distance
    alpha_ij.append(np.rad2deg(np.arctan2(xij[1], xij[2])))  # azimuth
    zij.append(np.rad2deg(np.arccos(xij[2] / sij)))
    # same operation for check
    sc = np.sqrt(pow(xij[0], 2) + pow(xij[1], 2) + pow(xij[2], 2))
    alpha_ij = np.rad2deg(np.arctan2(xij[1], xij[2]))
    zc = np.rad2deg(np.arccos(xij[2] / sij))
    print(sij)
    print(alpha_ij)
    print(zij)
    if (zc >= 0 and zc <= 90):
        cnt = cnt + 1
print(sij)
print(alpha_ij)
print(zij)
print(cnt)
