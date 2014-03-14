import numpy as np
from pylab import *

print 'reading file..'
d      = np.genfromtxt('worldcitiespop.txt',dtype='str',delimiter=',',skip_header=1)
cc     = d[:,0]   # country
city1  = d[:,1]   # name  
city2  = d[:,2]   # Name
popu   = d[:,4]   # population, convert to float later 
lat    = array(d[:,5],dtype=np.float32)
lon    = array(d[:,6],dtype=np.float32)

clist = ['nl','be','de','fr']

# filter - Quick and dirty manual loop
print 'filtering file..'
out = open('cities_eur.txt','w')
for i in range(len(lat)):
  if(i%float(int(len(lat)/10))==0):
    print i,len(lat)
  #if((lat[i]>40. and lat[i]<55) and (lon[i]>-5 and lon[i]<20)):
  if(cc[i] in clist):
    try:
      pop = float(popu[i])
      if(pop>200000):
        out.write('%s,%s,%.5f,%.5f\n'%(city1[i],city1[i][:3],lon[i],lat[i]))
    except:
      bla = 'aap'
out.close()

# test - read back in
#d      = np.genfromtxt('cities_eur.txt',dtype='str',delimiter=',')
#icao2  = d[:,0]
#lat2   = array(d[:,1],dtype=np.float32)
#lon2   = array(d[:,2],dtype=np.float32)
