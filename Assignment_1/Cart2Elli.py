import numpy as np

##parameters
#latitude - phi
LaD=59 #47
LaM=20 #0
LaS=59 #0
#longitude - lambda
LoD=18 #15
LoM=4 #0
LoS=10 #0

lat=LaD+LaM/60+LaS/3600
lon=LoD+LoM/60+LoS/3600
h=62+23
print(lat, lon ,h)
def Car2Ellip(lat, lon, h):
    a=6378137 # semi major axis
    fe=1 / 298.257222101
    b=a-fe*a # semi minor axis
    N=pow(a,2)/np.sqrt(pow(a,2)*pow(np.cos(np.deg2rad(lat)),2)+pow(b,2)*pow(np.sin(np.deg2rad(lat)),2))
    X=(N+h)*np.cos(np.deg2rad(lat))*np.cos(np.deg2rad(lon))
    Y=(N+h)*np.cos(np.deg2rad(lat))*np.sin(np.deg2rad(lon))
    Z=(pow(b/a,2)*N+h)*np.sin(np.deg2rad(lat))
    #print(X,Y,Z)
    return X,Y,Z

[X,Y,Z]=Car2Ellip(lat,lon,h)
#print(X,Y,Z)
##print (X,Y,Z)
#Part B - From global to local
#use the csv file
#X
#define another position to test the code 100m above the point

[X,Y,Z]=Car2Ellip(lat,lon,h)
#print(X,Y,Z)
#invert cart. back to ellipse

def Ellop2Car(X,Y,Z):
    a=6378137 # semi major axis
    fe=1 / 298.257222101
    b=a-fe*a # semi minor axis
    e2=(pow(a,2)-pow(b,2))/pow(a,2)
    p=np.sqrt(pow(X,2)+pow(Y,2))
    lat0=np.rad2deg(np.arctan2((Z/p)*(pow((1-e2),-1)),1))
    lat=0
    while(lat!=lat0):
        N0=pow(a,2)/np.sqrt(pow(a,2)*pow(np.cos(np.deg2rad(lat0)),2)+pow(b,2)*pow(np.sin(np.deg2rad(lat0)),2))
        h=p/np.cos(np.deg2rad(lat0))-N0
        lat=np.rad2deg(np.arctan2((Z/p)*pow(1-e2*(N0/(N0+h)),-1),1))
        lat0=lat
    lon=np.rad2deg(np.arctan2(Y,X))
    h=p/np.cos(np.deg2rad(lat0))-N0
    return lat0,lon,h


[lat,lon,h]=Ellop2Car(X,Y,Z)
print(lat,lon,h)
#start from global to local part