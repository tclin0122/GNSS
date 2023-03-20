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

a=6378137 # semi major axis
fe=1 / 298.257222101
b=a-fe*a # semi minor axis
#meter - Height above sea level
N=pow(a,2)/np.sqrt(pow(a,2)*pow(np.cos(np.deg2rad(lat)),2)+pow(b,2)*pow(np.sin(np.deg2rad(lat)),2))
print(N)
h=62+23
## converting calculation

X=(N+h)*np.cos(np.deg2rad(lat))*np.cos(np.deg2rad(lon))
Y=(N+h)*np.cos(np.deg2rad(lat))*np.sin(np.deg2rad(lon))
Z=(pow(b/a,2)*N+h)*np.sin(np.deg2rad(lat))
print (X,Y,Z)