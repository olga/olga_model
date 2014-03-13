# Prevents figs from showing on screen. 
# Call before pylab, etc.
import matplotlib
matplotlib.use('Agg')    

import numpy as np
from pylab import *
from mpl_toolkits.basemap import Basemap
import os
import gc
import sys
import urllib2
from copy import deepcopy
from scipy.ndimage.filters import gaussian_filter as gausf
import matplotlib.gridspec as gridspec
import matplotlib.image as image

# Custom modules
from readwrf import *
from colormaps import *
from readsounding import *
from sounding_v3 import *
from readopenair import *
from tools import *
from constants import *

# MOVE TO MATPLOTLIBRC
from matplotlib import rc
rc('font', size=9)
rc('legend', fontsize=8)

# -------------------------------------------------
# Function to create soundings
# -------------------------------------------------
def create_sounding(wrfout,domain,date,names,lons,lats,times):

  for name,lon,lat in zip(names,lons,lats):
    # Read WRF data
    d = readwrf_loc(wrfout,domain,lon,lat)
    sset = skewt_input()
    for t in times:
      for stype in range(2):
        print "sounding %s, type=%i, time=%i"%(name,stype,t)
        sset.stype  = stype
        sset.T      = d.T[t,:] 
        sset.Td     = d.Td[t,:]
        sset.p      = d.p[t,:] 
        sset.z      = d.zf[t,:]
        sset.u      = d.u[t,:]
        sset.v      = d.v[t,:]
        sset.name   = name
        sset.time   = d.datetime[t]
        sset.parcel = True
        sset.ps     = d.ps[t]
        sset.Ts     = d.T2[t]
        sset.rs     = d.q2[t]

        skewtlogp(sset)
      
        nameo = 'figures/'+ date + '/d' + str(domain) + '_sound' + str(stype) +'_' + name + '_' + str(t).zfill(2) + '.png'
        savefig(nameo)

        close('all')

# -------------------------------------------------
# Function to create (sort-of) meteograms
# -------------------------------------------------
def create_tser(wrfout,domain,date,names,lons,lats):

  levs = [0,1,2,3,4]
  wupd = cmap_discrete(wup,np.linspace(0,1,4))

  # Loop over requested locations
  for name,lon,lat in zip(names,lons,lats):
    # Read WRF data
    d = readwrf_loc(wrfout,domain,lon,lat)
    t = np.arange(0,24.001,1)

    fig = plt.figure(figsize=(6.5,6.0))
    #                   L    B    R    T    ws  hs
    fig.subplots_adjust(0.10,0.11,0.96,0.88,0.36,0.36)
    figtext(0.5,0.95,'%s [%.2fN, %.2fE]'%(name,lon,lat),size=9,ha='center')
    figtext(0.5,0.93,'%s'%(d.datetime[0]),size=8,ha='center')
    gs = gridspec.GridSpec(3,1,height_ratios=[4,1,1.5])

    ax=subplot(gs[0])
    ax.set_title('Updraft velocity and height',loc='left')
    zs = d.z[0,0]             # terrain height (lowest half level)
    wm = 3.5                  # scaling for colormap
    for i in range(t.size):
      bar(t[i]-0.35,d.zi[i],width=0.7,bottom=zs,color=wup((floor(d.wstar[i])+0.5)/wm),edgecolor='none')    
      bar(t[i]-0.4,d.ct[i]-d.zi[i],width=0.8,bottom=d.zi[i]+zs,color='k',alpha=0.3,edgecolor='none')    
    # Add sort-of colorbar
    wups = ([0.5,1.5,2.5,3.5])
    names = (['0-1 m/s','1-2 m/s','2-3 m/s','>3 m/s'])
    for wu,nam in zip(wups,names):
      scatter([-10],[300],color=wup(wu/wm),label=nam)
    legend(frameon=False,loc=2)  
    # finish things 
    xlim(0,24)
    ylim(0,3000)
    modplot(ax) 
    ylabel('z [m AMSL]')
    xticks(np.arange(0,24.001,2))

    ax=subplot(gs[1])
    ax.set_title('Cloud cover',loc='left')
    #ax.spines['left'].set_visible(False)
    #ax.get_yaxis().set_visible(False)
    pcolormesh(d.ccl,cmap=cm.bone_r,vmin=0,vmax=1)
    text(-0.4,0.5,'low',size=7,ha='right',va='center')
    text(-0.4,1.5,'middle',size=7,ha='right',va='center')
    text(-0.4,2.5,'high',size=7,ha='right',va='center')
    xlim(0,24)
    ylabel('z [m]')
    modplot(ax)
    ax.set_yticks([])
    xticks(np.arange(0,24.001,2))

    ax=subplot(gs[2])
    ax.set_title('Wind',loc='left')
    k500 = key_nearest(d.zf[0,:],500)
    k1000 = key_nearest(d.zf[0,:],1000)
    k2000 = key_nearest(d.zf[0,:],2000)
    barbs(t,0.5,d.u10*m2k,d.v10*m2k,length=5,linewidth=0.5,pivot='middle')  
    barbs(t,1.5,d.u[:,k500]*m2k,d.v[:,k500]*m2k,length=5,linewidth=0.5,pivot='middle')  
    barbs(t,2.5,d.u[:,k1000]*m2k,d.v[:,k1000]*m2k,length=5,linewidth=0.5,pivot='middle')  
    barbs(t,3.5,d.u[:,k2000]*m2k,d.v[:,k2000]*m2k,length=5,linewidth=0.5,pivot='middle')  
    text(-0.4,0.5,'10m',size=7,ha='right',va='center')
    text(-0.4,1.5,'500m',size=7,ha='right',va='center')
    text(-0.4,2.5,'1000m',size=7,ha='right',va='center')
    text(-0.4,3.5,'2000m',size=7,ha='right',va='center')
    xlim(0,24)
    modplot(ax)
    ax.set_yticks([])
    xticks(np.arange(0,24.001,2))
    xlabel('time UTC [h]')

    # Add logo :)
    img = image.imread('data/olga_lr.png')
    w=650;h=600
    figimage(img,10,6)
    figtext(0.08,0.013,'Open Limited-area Gliding Analysis. 6 x 6 km GFS-initiated WRF-ARW forecast [olga.vanstratum.com]',size=7,ha='left')


    name = 'figures/'+ date + '/d' + str(domain) + '_tser_' + name + '.png'
    savefig(name)
    savefig('test.png')


# -------------------------------------------------
# Function to create maps
# -------------------------------------------------
def create_maps(wrfout,domain,date,t0,t1,dt,variables,basemap,filter=False):

  print "creating maps"

  fsigma = 1.0  # std dev of gaussian filter size..

  # Read WRF output
  d  = readwrf_all(wrfout,domain=domain)

  makepfd = True

  # Get basemap coords cities
  if(domain==1):
    db     = np.genfromtxt('data/cities_eur_d1.txt',dtype='str',delimiter=',',skip_header=0)
  elif(domain==2):
    db     = np.genfromtxt('data/cities_eur_d2.txt',dtype='str',delimiter=',',skip_header=0)
  city   = db[:,0]
  citys  = db[:,1]  
  citlon = array(db[:,2],dtype=np.float32)  
  citlat = array(db[:,3],dtype=np.float32)  
  citlon,citlat = basemap(citlon,citlat)

  for t in range(t0,t1+1,dt):
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
        pfilt  = gausf(d.slps[t,:,:]/100.,3.0,mode='reflect')
        cf     = False #m.contourf(lon,lat,utot,levs2,extend='both',cmap=wnd)
        barbs  = m.barbs(lon[0::intv,0::intv],lat[0::intv,0::intv],\
                         d.U10[t,0::intv,0::intv],d.V10[t,0::intv,0::intv],\
                         length=5,barbcolor='r',flagcolor='r',linewidth=0.5)
        c      = m.contour(lon,lat,pfilt,levs1,linewidths=1.5,colors='k')
        cl     = clabel(c,levs1,fmt='%.1f')   
        doplot = True 
 
      # -------------------------------------------------
      # rain
      # -------------------------------------------------
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
      # Cloud cover
      # -------------------------------------------------
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

      if(var == 'clouds2'):
        title = 'Cloud cover -> NEEDS OPTIMIZATION'
        levs  = np.arange(0.00,1.,0.05)
        #intv  = 3
        cf    = m.contourf(lon,lat,d.calb[t,:,:]*6.,levs,alpha=1.,cmap=cloud)
        #lim   = 0.1
        doplot = True

 
      # -------------------------------------------------
      #   GLIDING SPECIFIC
      # -------------------------------------------------
 
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

      # -------------------------------------------------
      # TEST: z/L
      # -------------------------------------------------
      if(var == 'zol'):
        title = 'z/L [-]'
        levs  = np.arange(-100.,0.01,5.)
        data  = gausf((d.zi[t,:,:]/d.L[t,:,:]),fsigma,mode='reflect') \
                if filter else d.zi[t,:,:]/d.L[t,:,:]
        cf    = m.contourf(lon,lat,data,levs,extend='both',cmap=cent)
        doplot = True 
 
      # -------------------------------------------------
      # Top of updraft (dry convection)
      # -------------------------------------------------
      if(var == 'zidry'):
        title = 'Updraft height [m amsl]'
        levs  = np.arange(0,2500.01,200.)
        data  = gausf(d.zi[t,:,:]+d.hgt[:,:],fsigma,mode='reflect') \
                if filter else d.zi[t,:,:]+d.hgt[:,:]
        cf    = m.contourf(lon,lat,data,levs,extend='both',cmap=wup)
        doplot = True 
  
      # -------------------------------------------------
      # Cumulus depth
      # -------------------------------------------------
      if(var == 'cudepth'):
        title = 'Cumulus depth [m]'
        levs  = np.arange(0,2000,100)
        data  = gausf(d.zct[t,:,:]-d.zi[t,:,:],fsigma,mode='reflect') \
                if filter else d.zct[t,:,:]-d.zi[t,:,:]
        cf =  m.contourf(lon,lat,data,levs,extend='both',cmap=wupnl)
        doplot = True 

      # -------------------------------------------------
      # Potential flight distance
      # -------------------------------------------------
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
      # -------------------------------------------------
      if(doplot):

        # Plot terrain height contours
        if(False):
          levs = arange(50,4000.01,50.)
          contour(lon,lat,d.hgt,levs,cmap=cm.PuBu,alpha=0.5)

        # Plot cities
        if(True):
          for i in range(len(citys)):
            m.scatter(citlon[i],citlat[i],s=4,alpha=1.0)
            text(citlon[i],citlat[i],' '+citys[i],size=8,ha='left',va='center',color='0.1')

        # Plot airspace
        if(False):
          plotairspace(m)

        # If filled contour, add colorbar
        if(cf != False):
          pos = ax.get_position()
          l,b,w,h = pos.bounds
          cax = axes([l+w,b+0.1,0.02,h-0.2])
          cb=colorbar(cf,drawedges=False,cax=cax)
          cb.ax.tick_params(labelsize=8) 
          cb.outline.set_color('white')
        axes(ax)
       
        mc = '0.6' if var=="clouds2" else '0.2'
        cl = 1. if domain==1 else 1.5
        m.drawcoastlines(linewidth=cl,color=mc)
        m.drawcountries(linewidth=1,color=mc)
        m.drawmapboundary()
        
        if(domain==1):
          figtext(0.065,0.025,'OLGA: Open Limited-area Gliding Analysis. 18 x 18 km GFS-initiated WRF-ARW forecast [olga.vanstratum.com]',size=7,ha='left')
          #m.drawrivers(linewidth=0.5,color='#0066FF')
          #m.drawmeridians(arange(0, 360, 5))
          #m.drawparallels(arange(30, 60, 5))
        elif(domain==2):
          figtext(0.065,0.025,'OLGA: Open Limited-area Gliding Analysis. 6 x 6 km GFS-initiated WRF-ARW forecast [olga.vanstratum.com]',size=7,ha='left')
          m.drawrivers(linewidth=0.5,color='#0066FF')
          #m.drawmeridians(arange(0, 360, 5))
          #m.drawparallels(arange(30, 60, 5))
        subtitle = str(d.datetime[t]) + ' UTC'
        ax.set_title(title,loc='left')
        ax.set_title(subtitle,loc='right')

        # Add logo :)
        img = image.imread('data/olga_lr.png')
        w=650;h=600
        #figimage(img,w-50,h-45)
        figimage(img,w-45,6)
        #figimage(img,7,5)
 
        name = 'figures/'+ date + '/d' + str(domain) + '_' + var + '_' + str(t).zfill(3) + '.png'
        savefig(name)
    
      # Cleanup!
      fig.clf()
      plt.close()
      gc.collect()
 
  



