import numpy as np
from netCDF4 import Dataset
from pylab import *
from mpl_toolkits.basemap import Basemap

close('all')

class readgeom:
  def __init__(self,path):
    nc               = Dataset(path,'r')
    self.xlat        = nc.variables["XLAT_M"][0,:,:]
    self.xlon        = nc.variables["XLONG_M"][0,:,:]
    self.hgt_m       = nc.variables["HGT_M"][0,:,:]

def plotdomain(mapin,lats,lons):
  lat1 = lats[0,0]
  lat2 = lats[-1,0]
  lon1 = lons[0,0]
  lon2 = lons[0,-1]
  lon3 = lons[-1,0]
  lon4 = lons[-1,-1]

  x1,y1=mapin(lon1,lat1)
  x2,y2=mapin(lon2,lat1)
  x3,y3=mapin(lon3,lat2)
  x4,y4=mapin(lon4,lat2)

  mapin.plot([x1,x2],[y1,y1],'w-',linewidth=2)
  mapin.plot([x2,x4],[y1,y3],'w-',linewidth=2)
  mapin.plot([x4,x3],[y3,y3],'w-',linewidth=2)
  mapin.plot([x3,x1],[y3,y1],'w-',linewidth=2)

h  = readgeom('back.nc')
d1 = readgeom('../WPS/geo_em.d01.nc')
d2 = readgeom('../WPS/geo_em.d02.nc')

fig = plt.figure(figsize=(9.4,9))
#                   L    B    R    T    ws  hs
fig.subplots_adjust(0.08,0.05,.99,0.95,0.2,0.08)
m = Basemap(width=2500000,height=2500000,
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=10.,projection='lcc',\
            lat_1=52.,lat_2=52.,lat_0=52.,lon_0=5.5)

ax = fig.add_axes([0.00,0.05,0.95,0.85])
x,y = m(h.xlon[:,:],h.xlat[:,:])
cf = m.pcolormesh(x,y,h.hgt_m[:,:],cmap=cm.RdBu_r)

plotdomain(m,d1.xlat,d1.xlon)
plotdomain(m,d2.xlat,d2.xlon)
 
pos = ax.get_position()
l,b,w,h = pos.bounds
cax = axes([l+w-0.05,b+0.1,0.02,h-0.2])
colorbar(cf,drawedges=False,cax=cax)
axes(ax)

m.drawcoastlines(linewidth=1,color='0.')
m.drawcountries(linewidth=1,color='0.')
m.drawmapboundary()
m.drawmeridians(arange(0, 360, 10))
m.drawparallels(arange(30, 60, 10))


