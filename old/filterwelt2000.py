import numpy as np
from pylab import *

d      = np.genfromtxt('welt2000.txt',dtype='str',delimiter=',',skip_header=22)
name   = d[:,1]
ptype  = d[:,2]
lat    = array(d[:,3],dtype=np.float32)
lon    = array(d[:,4],dtype=np.float32)
icao   = d[:,15]

# filter - Quick and dirty manual loop
out = open('welt2000f.txt','w')
for i in range(len(lat)):
  if(lat[i]>40. and lat[i]<55 and lon[i]>-5 and lon[i]<20 and ptype[i]=='APT'):
    out.write('%s,%.5f,%.5f\n'%(icao[i],lon[i],lat[i]))
out.close()

# test - read back in
d      = np.genfromtxt('welt2000f.txt',dtype='str',delimiter=',')
icao2  = d[:,0]
lat2   = array(d[:,1],dtype=np.float32)
lon2   = array(d[:,2],dtype=np.float32)
