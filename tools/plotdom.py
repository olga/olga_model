import numpy as np
from pylab import *
from math import radians, cos, sin, asin, sqrt
from mpl_toolkits.basemap import Basemap

# ====== SETTINGS =======
ndomains = 2
#WPSPath = '../domain_test/WPS' 
WPSPath = '../OLGA_NL/WPS' 
# =======================

def get_nc_obj(file):
    # Check if the netcdf4-python module is available:
    try:
        from netCDF4 import Dataset
        wrfin = Dataset(file,'r')
        netcdf4py = True 
    except:
        netcdf4py = False 

    # Fallback option: NetCDF from Scientific.IO
    try:
        from Scientific.IO import NetCDF
        wrfin = NetCDF.NetCDFFile(file,'r')
        netcdfsci = True
    except:
        netcdfsci = False

    if(netcdf4py==False and netcdfsci==False):
        sys.exit('No NetCDF module available!')
    else:
        return wrfin

def haversine(lon1, lat1, lon2, lat2):
    lon1,lat1,lon2,lat2 = map(radians,[lon1,lat1,lon2,lat2])
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    return 6371 * c

class readgeom:
    def __init__(self,path):
        print('reading %s'%path)
        nc               = get_nc_obj(path)
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

    mapin.plot([x1,x2],[y1,y1],'r-',linewidth=2)
    mapin.plot([x2,x4],[y1,y3],'r-',linewidth=2)
    mapin.plot([x4,x3],[y3,y3],'r-',linewidth=2)
    mapin.plot([x3,x1],[y3,y1],'r-',linewidth=2)

# Read the geo_em domain files
domains = []
for n in range(ndomains):
    domains.append(readgeom('%s/geo_em.d%02i.nc'%(WPSPath,n+1)))

# Get some crude info of domain locations
centlat = np.average(domains[0].xlat)
centlon = np.average(domains[0].xlon)
width   = 1.05* haversine(domains[0].xlon.min(),centlat,domains[0].xlon.max(),centlat) * 1000 
height  = 1.05* haversine(centlon, domains[0].xlat.min(), centlon, domains[0].xlat.max()) * 1000

# Plot maps
if(True):
    close('all')

    fig = plt.figure(figsize=(10,0.85*10/width*height))
    #                   L    B    R    T    ws  hs
    fig.subplots_adjust(0.08,0.05,.99,0.95,0.2,0.08)
    ax = subplot(111)
    m = Basemap(width=width,height=height,
                rsphere=(6378137.00,6356752.3142),\
                resolution='c',area_thresh=10.,projection='lcc',\
                lat_1=centlat,lat_2=centlat,lat_0=centlat,lon_0=centlon)
   
    # Plot outline domain(s)
    for n in range(ndomains):
        plotdomain(m,domains[n].xlat, domains[n].xlon)
    
        # Plot topography of outer domain
        if(n==0):
            x,y = m(domains[n].xlon[:,:],domains[n].xlat[:,:])
            cf = m.pcolormesh(x,y,domains[n].hgt_m[:,:],cmap=cm.RdBu_r)
            cf = m.pcolormesh(x,y,domains[n].hgt_m[:,:],cmap=cm.RdBu_r)
            colorbar()
    
    m.drawcoastlines(linewidth=1,color='0.')
    m.drawcountries(linewidth=1,color='0.')
    m.drawmapboundary()
    m.drawmeridians(np.arange(0, 360, 5),labels=[True,False,False,True])
    m.drawparallels(np.arange(30, 60, 5),labels=[True,False,False,True])

    savefig('domains.png')
    #savefig('domains.pdf')
