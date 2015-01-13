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
from tools import *

## Function to create maps
def create_maps(olga,wrfout,dom,times):
    filter = True
    fsigma = 0.5     # std dev of gaussian filter size..

    # Read WRF output
    d  = readwrf_all(olga,wrfout,times[0],times[-1]) 

    # If any map parameters are set to -1, try to determine best setting
    if(olga.map_lat[dom] == -1):
        map_lat = np.average(d.lat)
    else:
        map_lat = olga.map_lat[dom]

    if(olga.map_lon[dom] == -1):
        map_lon = np.average(d.lon)
    else:
        map_lon = olga.map_lat[dom]

    if(olga.map_width[dom] == -1):
        map_width =  1.05* haversine(d.lon.min(), map_lat, d.lon.max(), map_lat) * 1000
    else:
        map_width = olga.map_width[dom]

    if(olga.map_height[dom] == -1):
        map_height =  1.05* haversine(map_lon, d.lat.min(), map_lon, d.lat.max()) * 1000
    else:
        map_height = olga.map_height[dom]

    # Setup basemap only once and (deep)copy later for performance
    basem = Basemap(width=map_width,height=map_height,
                    rsphere=(6378137.00,6356752.3142),\
                    resolution=olga.map_res[dom],area_thresh=10.,projection='lcc',\
                    lat_1=map_lat,lat_0=map_lat,lon_0=map_lon)

    olga_logo = pl.matplotlib.image.imread(olga.olgaRoot+'include/olga_lr.png')

    # Get basemap coords cities
    #if(domain==1):
    #    db = np.genfromtxt(olgaRoot+'include/cities_eur_d1.txt',dtype='str',delimiter=',',skip_header=0)
    #elif(domain==2):
    #    db = np.genfromtxt(olgaRoot+'include/cities_eur_d2.txt',dtype='str',delimiter=',',skip_header=0)
    #city   = db[:,0]
    #citys  = db[:,1]  
    #citlon = np.array(db[:,2],dtype=np.float32)  
    #citlat = np.array(db[:,3],dtype=np.float32)  
    #citlon,citlat = basem(citlon,citlat)

    # Loop over differnt time steps 
    for t in times:
        for var in olga.map_vars[dom]:
            #print 'time = %s, var = %s'%(d.datetime[t],var)
  
            fig        = pl.figure(figsize=(6.5,6.0))
            m          = deepcopy(basem)
            axloc      = [0.03,0.05,0.85,0.87]  # left,bottom,width,height
            ax         = fig.add_axes(axloc)
            lon,lat    = m(d.lon,d.lat)
            doplot     = False

            # -------------------------------------------------
            # GENERAL  METEOROLOGY
            # -------------------------------------------------

            # -------------------------------------------------
            # sea level pressure and 10m wind
            # -------------------------------------------------
            if(var == 'slpwind'):
                title  = 'Sea level pressure (hPa) + 10m wind (kts)'
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
            # stream lines of wind
            # -------------------------------------------------
            if(var[:4] == 'wind'):
                if(var=='wind10m'):
                    title  = '10m wind (kts)'
                    ufilt  = gausf(d.U10[t,:,:],fsigma,mode='reflect')
                    vfilt  = gausf(d.V10[t,:,:],fsigma,mode='reflect')
                    scalefac = 25.
                elif(var=='wind1000m'):
                    title  = '1000m wind'
                    ufilt  = gausf(d.U1000[t,:,:],fsigma,mode='reflect')
                    vfilt  = gausf(d.V1000[t,:,:],fsigma,mode='reflect')
                    scalefac = 35.

                utot   = ((ufilt[:,:]*m2k)**2.+(vfilt[:,:]*m2k)**2.)**0.5
                umax   = 40. #np.ceil((((d.U10[:,:,:]*m2k)**2.+(d.V10[:,:,:]*m2k)**2.)**0.5).max()*10.)/10.
                levs2  = np.arange(0,umax+0.01,umax/40.)  # Wind levels
                cf     = m.contourf(lon,lat,utot,levs2,extend='both',cmap=wnd)
                units  = 'kts'
                stream = m.streamplot(lon[0,:],lat[:,0],\
                                 ufilt[:,:],vfilt[:,:],density=3,linewidth=utot/scalefac,color='k')
                doplot = True 
 
            # -------------------------------------------------
            # rain
            # -------------------------------------------------
            if(var == 'rr'):
                title = 'Precipitation over last hour (mm, o=convective)'
                if(t==0):
                    rr  = d.rr_mp[t,:,:]+d.rr_con[t,:,:]
                    rrc = d.rr_con[t,:,:]
                else:
                    rr  = (d.rr_mp[t,:,:]+d.rr_con[t,:,:])-(d.rr_mp[t-1,:,:]+d.rr_con[t-1,:,:])  
                    rrc = d.rr_con[t,:,:]-d.rr_con[t-1,:,:]  
                levs  = np.arange(0.1,30.0001,0.5)
                cf    = m.contourf(lon,lat,rr,levs,alpha=1.,extend='both',cmap=rain3)
                units = 'mm'
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
                title = 'Cloud cover fraction (o=low  -=mid  |=high)'
                levs  = np.arange(0.00,1.001,0.05)
                intv  = 3
                cf    = m.contourf(lon,lat,d.ccsum[t,:,:],levs,alpha=1.,cmap=cld)
                units = 'fraction'
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

            # -------------------------------------------------
            # Shortwave downwelling radiation
            # TO-DO: make it potential (i.e. max potential swd / actual swd)
            # -------------------------------------------------
            if(var == 'swd'):
                title = 'Shortwave incoming radiation (W/m2)'
                levs  = np.arange(0,1000.01,50)
                data  = gausf(d.swd[t,:,:],fsigma,mode='reflect') \
                        if filter else d.swd[t,:,:]
                cf    = m.contourf(lon,lat,data,levs,extend='both',cmap=wup)
                units = 'W/m2'
                doplot = True 
 
            # -------------------------------------------------
            # GLIDING SPECIFIC
            # -------------------------------------------------

            # -------------------------------------------------
            # Updraft velocity wstar
            # -------------------------------------------------
            if(var == 'wstar'):
                title = 'Updraft velocity (m/s)'
                levs  = np.arange(0,3.001,0.5)
                data  = gausf(d.wstar[t,:,:],fsigma,mode='reflect') \
                        if filter else d.wstar[t,:,:]
                cf    = m.contourf(lon,lat,data,levs,extend='both',cmap=wup)
                units = 'm/s'
                doplot = True 

            # -------------------------------------------------
            # TEST: z/L as indicator for cloud streets / roll convection
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
                title = 'Updraft height [m AMSL]'
                levs  = np.arange(0,2500.01,300.)
                data  = gausf(d.zi[t,:,:]+d.hgt[:,:],fsigma,mode='reflect') \
                        if filter else d.zi[t,:,:]+d.hgt[:,:]
                cf    = m.contourf(lon,lat,data,levs,extend='both',cmap=wup)
                units = 'm AMSL'
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
                units = 'm'
                doplot = True 

            # -------------------------------------------------
            # Potential flight distance per day
            # -------------------------------------------------
            if((var == 'pfd') and t < np.size(d.date_PFD)):
                title = 'Potential flight distance [km]'
                levs  = np.arange(0,900.1,100)
                data  = gausf(d.PFD[t,:,:],fsigma,mode='reflect') \
                        if filter else d.PFD[t,:,:]
                cf    =  m.contourf(lon,lat,data,levs,extend='both',cmap=wup)
                units = 'km'
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
                #if(True):
                #    for i in range(len(citys)):
                #        m.scatter(citlon[i],citlat[i],s=4,alpha=1.0)
                #        pl.text(citlon[i],citlat[i],' '+citys[i],size=8,ha='left',va='center',color='0.1')

                # Plot airspace
                #if(False):
                #    plotairspace(m)

                # If filled contour, add colorbar
                if(cf != False):
                    pos = ax.get_position()
                    l,b,w,h = pos.bounds
                    cax = pl.axes([l+w,b+0.1,0.02,h-0.2])
                    cb=pl.colorbar(cf,drawedges=False,cax=cax)
                    cb.ax.tick_params(labelsize=8) 
                    cb.outline.set_color('white')
                    cb.set_label(units)
                pl.axes(ax)
           
                m.drawcoastlines(linewidth=1,color='0.5')
                m.drawcountries(linewidth=1,color='0.5')
                #m.drawrivers(linewidth=0.5,color='#0066FF')
                m.drawmapboundary()
               
                # Add description variables, units, OLGA label, logo, etc. 
                pl.figtext(0.065,0.025,'OLGA: %s GFS-initiated WRF-ARW forecast [olga.vanstratum.com]'\
                           %olga.map_desc[dom],size=6,ha='left')
                if(var=='pfd'):
                    subtitle = 'cumulative distance over: %s'%str(d.date_PFD[t]) 
                else:
                    subtitle = 'valid: %s UTC'%str(d.datetime[t])

                # Add OLGA logo
                ax.set_title(title,loc='left')
                ax.set_title(subtitle,loc='right')
                w=650 ; h=600
                pl.figimage(olga_logo,w-45,6)
 
                if(var=='pfd'):
                    tmp    = '%04i'%(t*24.) 
                else:
                    xtime  = d.time[t] / 3600.
                    hour   = int(np.floor(xtime))
                    minute = int((xtime - hour) * 60.)
                    tmp    = str(hour).zfill(3) + str(minute).zfill(2)

                name = '%s%04i%02i%02i_t%02iz/%s_d%i_%s.png'%(olga.figRoot,olga.year,olga.month,olga.day,olga.cycle,tmp,dom+1,var)
                pl.savefig(name)
        
            # Cleanup!
            fig.clf()
            pl.close()
            gc.collect()

