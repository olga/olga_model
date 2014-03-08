import numpy as np
from pylab import *
from mpl_toolkits.basemap import Basemap
import os
import gc
import sys
import urllib2
from copy import deepcopy
from scipy.ndimage.filters import gaussian_filter as gausf

# Custom modules
from readwrf import *
from colormaps import *
from readsounding import *
from soundingv2 import *
from readopenair import *

# MOVE TO MATPLOTLIBRC
from matplotlib import rc
rc('font', size=9)
rc('legend', fontsize=12)

# -------------------------------------------------
#   Function to create maps
# -------------------------------------------------
def create_maps(wrfout,domain,date,t0,t1,variables,basemap,filter=False):

  fsigma = 0.5  # std dev of gaussian filter size..

  # Read WRF output
  d  = readwrf_all(wrfout,domain=domain)

  makepfd = True

  # Get basemap coords cities
  if(domain==1):
    db     = np.genfromtxt('cities_eur_d1.txt',dtype='str',delimiter=',',skip_header=0)
  elif(domain==2):
    db     = np.genfromtxt('cities_eur_d2.txt',dtype='str',delimiter=',',skip_header=0)
  city   = db[:,0]
  citys  = db[:,1]  
  citlon = array(db[:,2],dtype=np.float32)  
  citlat = array(db[:,3],dtype=np.float32)  
  citlon,citlat = basemap(citlon,citlat)

  for t in range(t0,t1+1,1):
    for var in variables:
      print 'time = %s, var = %s'%(d.datetime[t],var)
   
      fig     = plt.figure(figsize=(6.5,6.0))
      m       = deepcopy(basemap)
      axloc   = [0.03,0.05,0.85,0.87]  # left,bottom,width,height
      ax      = fig.add_axes(axloc)
      lon,lat = m(d.lon,d.lat)
      doplot  = False

      # -------------------------------------------------
      #   GENERAL  METEOROLOGY
      # -------------------------------------------------

      # sea level pressure and 10m wind
      # -------------------------------------------------
      if(var == 'slpwind'):
        title  = 'Sea level pressure [hPa] + 10m wind [kts]'
        utot   = (d.U10[t,:,:]**2.+d.V10[t,:,:]**2.)**0.5
        umax   = np.ceil(((d.U10[:,:,:]**2.+d.V10[:,:,:]**2.)**0.5).max()*10.)/10.
        intv   = 4    # interval of wind barbs
        levs1  = np.arange(800,1200,2.5)     # SLP levels
        levs2  = np.arange(0,umax,umax/50.)  # Wind levels
        pfilt  = gaussian_filter(d.slps[t,:,:]/100.,3.0,mode='reflect')
        cf     = False #m.contourf(lon,lat,utot,levs2,extend='both',cmap=wnd)
        barbs  = m.barbs(lon[0::intv,0::intv],lat[0::intv,0::intv],\
                         d.U10[t,0::intv,0::intv],d.V10[t,0::intv,0::intv],\
                         length=5,barbcolor='r',flagcolor='r',linewidth=0.5)
        c      = m.contour(lon,lat,pfilt,levs1,linewidths=1.5,colors='k')
        cl     = clabel(c,levs1,fmt='%.1f')   
        doplot = True 
 
      # -------------------------------------------------
      # rain
      if(var == 'rr'):
        title = 'Precipitation over last hour [mm] (o=convective)'
        if(t==0):
          rr  = d.rr_mp[t,:,:]+d.rr_con[t,:,:]
          rrc = d.rr_con[t,:,:]
        else:
          rr  = (d.rr_mp[t,:,:]+d.rr_con[t,:,:])-(d.rr_mp[t-1,:,:]+d.rr_con[t-1,:,:])  
          rrc = d.rr_con[t,:,:]-d.rr_con[t-1,:,:]  
        levs  = arange(0.1,30.0001,0.5)
        cf    = m.contourf(lon,lat,rr,levs,alpha=1.,extend='both',cmap=rain3)
        lim   = 0.1 ; intv  = 3
        for i in range(0,d.nlon,intv):
          for j in range(0,d.nlat,intv):
            if(rrc[i,j] > lim):
              text(lon[i,j],lat[i,j],'o',size=6,ha='center',va='center',color='0.1')
        doplot = True

      # -------------------------------------------------
      # Clouds
      if(var == 'clouds'):
        title = 'Cloud cover (o=low  -=mid  |=high)'
        levs  = np.arange(0.00,1.001,0.05)
        intv  = 3
        cf    = m.contourf(lon,lat,d.ccsum[t,:,:],levs,alpha=1.,cmap=cld)
        lim   = 0.1
        for i in range(0,d.nlon,intv):
          for j in range(0,d.nlat,intv):
            if(d.cclow[t,i,j] > lim):
              text(lon[i,j],lat[i,j],'o',size=8,ha='center',va='center',color='0.1')
            if(d.ccmid[t,i,j] > lim):
              text(lon[i,j],lat[i,j],'-',size=9,ha='center',va='center',color='0.2')
            if(d.cchig[t,i,j] > lim):
              text(lon[i,j],lat[i,j],'|',size=8,ha='center',va='center',color='0.1')
        doplot = True
 
      # -------------------------------------------------
      #   GLIDING SPECIFIC
      # -------------------------------------------------
 
      # Updraft velocity wstar
      # -------------------------------------------------
      if(var == 'wstar'):
        title = 'Updraft velocity (w*) [m/s]'
        levs  = np.arange(0,3.001,0.5)
        data  = gausf(d.wstar[t,:,:],fsigma,mode='reflect') \
                if filter else d.wstar[t,:,:]
        cf    = m.contourf(lon,lat,data,levs,extend='both',cmap=wup)
        doplot = True 
 
      # Top of updraft (dry convection)
      # -------------------------------------------------
      if(var == 'zidry'):
        title = 'Updraft height [m]'
        levs  = np.arange(0,2500.01,150)
        data  = gausf(d.zi[t,:,:]+d.hgt[:,:],fsigma,mode='reflect') \
                if filter else d.zi[t,:,:]+d.hgt[:,:]
        cf    = m.contourf(lon,lat,data,levs,extend='both',cmap=wup)
        doplot = True 
  
      # -------------------------------------------------
      # Cumulus depth
      if(var == 'cudepth'):
        title = 'Cumulus depth [m]'
        levs  = np.arange(0,2000,100)
        data  = gausf(d.zct[t,:,:]-d.zi[t,:,:],fsigma,mode='reflect') \
                if filter else d.zct[t,:,:]-d.zi[t,:,:]
        cf =  m.contourf(lon,lat,data,levs,extend='both',cmap=wupnl)
        doplot = True 

      # -------------------------------------------------
      # Potential flight distance
      if(var == 'pfd' and makepfd):
        pfd   = getpfd(wrfout)
        title = 'Potential flight distance [km]'
        levs  = np.arange(0,1000.1,100)
        data  = gausf(pfd.pfd[-1,:,:],fsigma,mode='reflect') \
                if filter else pfd.pfd[-1,:,:]
        cf    =  m.contourf(lon,lat,data,levs,extend='both',cmap=wup)
        makepfd = False
        doplot = True 

      # -------------------------------------------------
      # Finish plot!
      if(doplot):
        if(False):
          levs = arange(50,4000.01,50.)
          contour(lon,lat,d.hgt,levs,cmap=cm.PuBu,alpha=0.5)

        if(True):
          for i in range(len(citys)):
            m.scatter(citlon[i],citlat[i],s=4,alpha=1.0)
            text(citlon[i],citlat[i],' '+citys[i],size=8,ha='left',va='center',color='0.1')

        if(True):
          plotairspace(m)

        if(cf != False):
          pos = ax.get_position()
          l,b,w,h = pos.bounds
          cax = axes([l+w,b+0.1,0.02,h-0.2])
          cb=colorbar(cf,drawedges=False,cax=cax)
          cb.ax.tick_params(labelsize=8) 
          cb.outline.set_color('white')
        axes(ax)
        
        m.drawcoastlines(linewidth=1.5,color='0.2')
        m.drawcountries(linewidth=1,color='0.2')
        m.drawmapboundary()
        
        if(domain==1):
          figtext(0.065,0.025,'18 x 18 km GFS-initiated WRF-ARW forecast [www.dummy.org]',size=7,ha='left')
          m.drawmeridians(arange(0, 360, 5))
          m.drawparallels(arange(30, 60, 5))
        elif(domain==2):
          figtext(0.065,0.025,'6 x 6 km GFS-initiated WRF-ARW forecast [www.dummy.org]',size=7,ha='left')
          m.drawrivers(linewidth=0.5,color='#0066FF')
          m.drawmeridians(arange(0, 360, 5))
          m.drawparallels(arange(30, 60, 5))
        subtitle = str(d.datetime[t]) + ' UTC'
        ax.set_title(title,loc='left')
        ax.set_title(subtitle,loc='right')
 
        name = 'figures/'+ date + '/d' + str(domain) + '_' + var + '_' + str(t).zfill(3) + '.png'
        savefig(name)
    
      # Cleanup!
      fig.clf()
      plt.close()
      gc.collect()
 
  



