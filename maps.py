#
# Copyright (c) 2013-2014 Bart van Stratum (bart@vanstratum.com)
# 
# This file is part of OLGA.
# 
# OLGA is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# OLGA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with OLGA.  If not, see <http://www.gnu.org/licenses/>.
#

import gc
import numpy as np
import pylab as pl
from mpl_toolkits.basemap import Basemap
from scipy.ndimage.filters import gaussian_filter as gausf

from readwrf import *
from colormaps import *

## Function to create maps
def create_maps(wrfout,domain,date,t0,t1,dt,variables,filter=False):

    # Setup basemap only once and (deep)copy later for performance
    if(domain==1):
        basem = Basemap(width=1800000,height=1800000,
                        rsphere=(6378137.00,6356752.3142),\
                        resolution='l',area_thresh=10.,projection='lcc',\
                        lat_1=51.4,lat_2=51.4,lat_0=51.4,lon_0=5.5)
    elif(domain==2):
        basem = Basemap(width=590000,height=590000,
                        rsphere=(6378137.00,6356752.3142),\
                        resolution='i',area_thresh=10.,projection='lcc',\
                        lat_1=51.3,lat_2=51.3,lat_0=51.3,lon_0=6.7)

    fsigma = 1.0     # std dev of gaussian filter size..
    iPFD   = 0       # index of pfd

    # Read WRF output
    d  = readwrf_all(wrfout,domain=domain)  

    # Get basemap coords cities
    if(domain==1):
        db = np.genfromtxt('data/cities_eur_d1.txt',dtype='str',delimiter=',',skip_header=0)
    elif(domain==2):
        db = np.genfromtxt('data/cities_eur_d2.txt',dtype='str',delimiter=',',skip_header=0)
    city   = db[:,0]
    citys  = db[:,1]  
    citlon = np.array(db[:,2],dtype=np.float32)  
    citlat = np.array(db[:,3],dtype=np.float32)  
    citlon,citlat = basem(citlon,citlat)

    # Loop over differnt time steps 
    for t in range(t0,t1+1,dt):
        for var in variables:
            print 'time = %s, var = %s'%(d.datetime[t],var)
  
            fig     = pl.figure(figsize=(6.5,6.0))
            m       = deepcopy(basem)
            axloc   = [0.03,0.05,0.85,0.87]  # left,bottom,width,height
            ax      = fig.add_axes(axloc)
            lon,lat = m(d.lon,d.lat)
            doplot  = False

            # GENERAL  METEOROLOGY
            # -------------------------------------------------
            # sea level pressure and 10m wind
            # -------------------------------------------------
            if(var == 'slpwind'):
                title  = 'Sea level pressure [hPa] + 10m wind [kts]'
                utot   = ((d.U10[t,:,:]*m2k)**2.+(d.V10[t,:,:]*m2k)**2.)**0.5
                umax   = np.ceil((((d.U10[:,:,:]*m2k)**2.+(d.V10[:,:,:]*m2k)**2.)**0.5).max()*10.)/10.
                intv   = 4 # interval of wind barbs
                levs1  = np.arange(800,1200,2.5)     # SLP levels
                levs2  = np.arange(0,umax,umax/50.)  # Wind levels
                pfilt  = gausf(d.slps[t,:,:]/100.,3.0,mode='reflect')
                cf     = False #m.contourf(lon,lat,utot,levs2,extend='both',cmap=wnd)
                barbs  = m.barbs(lon[0::intv,0::intv],lat[0::intv,0::intv],\
                                 d.U10[t,0::intv,0::intv]*m2k,d.V10[t,0::intv,0::intv]*m2k,\
                                 length=5,barbcolor='r',flagcolor='r',linewidth=0.5)
                c      = m.contour(lon,lat,pfilt,levs1,linewidths=1.5,colors='k')
                cl     = pl.clabel(c,levs1,fmt='%.1f')   
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
                levs  = np.arange(0.1,30.0001,0.5)
                cf    = m.contourf(lon,lat,rr,levs,alpha=1.,extend='both',cmap=rain3)
                lim   = 0.1 ; intv  = 3
                for i in range(0,d.nlon,intv):
                    for j in range(0,d.nlat,intv):
                        if(rrc[i,j] > lim):
                            pl.text(lon[i,j],lat[i,j],'o',size=6,ha='center',va='center',color='0.1')
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
                            pl.text(lon[i,j],lat[i,j],'o',size=8,ha='center',va='center',color='0.1')
                        if(d.ccmid[t,i,j] > lim):
                            pl.text(lon[i,j],lat[i,j],'-',size=9,ha='center',va='center',color='0.2')
                        if(d.cchig[t,i,j] > lim):
                            pl.text(lon[i,j],lat[i,j],'|',size=8,ha='center',va='center',color='0.1')
                doplot = True

            if(var == 'clouds2'):
                title = 'Cloud cover -> NEEDS OPTIMIZATION'
                levs  = np.arange(0.00,1.,0.05)
                cf    = m.contourf(lon,lat,d.calb[t,:,:]*6.,levs,alpha=1.,cmap=cloud)
                doplot = True

 
            # GLIDING SPECIFIC
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
                levs  = np.arange(0,2500.01,300.)
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
            #if(var == 'pfd' and t*dt%24==0 and t!=t1):
            #    pfd   = getpfd(wrfout,t,t+(24/dt))
            #    title = 'Potential flight distance [km]'
            #    levs  = np.arange(0,900.1,100)
            #    data  = gausf(pfd.pfd[:,:],fsigma,mode='reflect') \
            #            if filter else pfd.pfd[:,:]
            #    cf    =  m.contourf(lon,lat,data,levs,extend='both',cmap=wup)
            #    makepfd = False
            #    doplot = True 

            if((var == 'pfd') and (t in d.tPFD)):
                title = 'Potential flight distance [km]'
                levs  = np.arange(0,900.1,100)
                data  = gausf(d.PFD[iPFD,:,:],fsigma,mode='reflect') \
                        if filter else d.PFD[iPFD,:,:]
                cf    =  m.contourf(lon,lat,data,levs,extend='both',cmap=wup)
                iPFD  += 1
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
                        pl.text(citlon[i],citlat[i],' '+citys[i],size=8,ha='left',va='center',color='0.1')

                # Plot airspace
                if(False):
                    plotairspace(m)

                # If filled contour, add colorbar
                if(cf != False):
                    pos = ax.get_position()
                    l,b,w,h = pos.bounds
                    cax = pl.axes([l+w,b+0.1,0.02,h-0.2])
                    cb=pl.colorbar(cf,drawedges=False,cax=cax)
                    cb.ax.tick_params(labelsize=8) 
                    cb.outline.set_color('white')
                pl.axes(ax)
             
                mc = '0.6' if var=="clouds2" else '0.2'
                cl = 1. if domain==1 else 1.5
                m.drawcoastlines(linewidth=cl,color=mc)
                m.drawcountries(linewidth=1,color=mc)
                m.drawmapboundary()
                
                if(domain==1):
                    pl.figtext(0.065,0.025,'OLGA: Open Limited-area Gliding Analysis. 18 x 18 km GFS-initiated WRF-ARW forecast [olga.vanstratum.com]',size=7,ha='left')
                elif(domain==2):
                    pl.figtext(0.065,0.025,'OLGA: Open Limited-area Gliding Analysis. 6 x 6 km GFS-initiated WRF-ARW forecast [olga.vanstratum.com]',size=7,ha='left')
                    m.drawrivers(linewidth=0.5,color='#0066FF')
                subtitle = str(d.datetime[t]) + ' UTC'
                ax.set_title(title,loc='left')
                ax.set_title(subtitle,loc='right')

                # Add logo :)
                img = pl.matplotlib.image.imread('data/olga_lr.png')
                w=650 ; h=600
                pl.figimage(img,w-45,6)
 
                tmp = '%02i_%02i'%(np.floor(t0+t*d.dt),(t0+t*d.dt-np.floor(t0+t*d.dt))*30.)
                name = 'figures/'+ date + '/' + tmp + '_d' + str(domain) + '_' + var + '_' + '.png'
                pl.savefig(name)
        
            # Cleanup!
            fig.clf()
            pl.close()
            gc.collect()

