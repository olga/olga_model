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
from matplotlib.colors import ListedColormap
from scipy.ndimage.filters import gaussian_filter as gausf

# Custom modules
from readwrf import *
from colormaps import *
from tools import *

"""
* Global settings: make property of settings?
"""
fig_width_px  = 650  # width of figure [px]
fig_dpi       = 100  # resolution in [px/in]
mar_bottom_px = 20   # bottom margin [px]
mar_top_px    = 25   # top margin in [px]
mar_left_px   = 12   # left margin in [px]
mar_right_px  = 80   # right margin in [px]

""" 
* Set the map properties, either from
* settings, or estimate from WRF output 
""" 
def getDomain(olga, dom, wrf):
    if(olga.map_lat[dom] == -1):
        map_lat = np.average(wrf.lat)
    else:
        map_lat = olga.map_lat[dom]

    if(olga.map_lon[dom] == -1):
        map_lon = np.average(wrf.lon)
    else:
        map_lon = olga.map_lon[dom]

    if(olga.map_width[dom] == -1):
        map_width =  1.05 * haversine(wrf.lon.min(), map_lat, wrf.lon.max(), map_lat) * 1000
    else:
        map_width = olga.map_width[dom]

    if(olga.map_height[dom] == -1):
        map_height =  1.05 * haversine(map_lon, wrf.lat.min(), map_lon, wrf.lat.max()) * 1000
    else:
        map_height = olga.map_height[dom]

    return (map_lat, map_lon, map_width, map_height)

""" 
* Calculate the figure properties to ensure:
* 1. Prescribed figure width in pixels
* 2. Fixed margins for labels, colorbar, etc 
"""
def getFigureProperties(map_height, map_width):
    # Aspect ratio of map
    map_aspect    = (float(map_height) / float(map_width))
   
    # Absolute size map
    map_width_px  = fig_width_px - mar_left_px - mar_right_px # [px]
    map_height_px = map_width_px * map_aspect # [px]

    # Relative left spacing and map width
    ax_left       = float(mar_left_px)  / float(fig_width_px) 
    ax_width      = float(map_width_px) / float(fig_width_px)
   
    # Figure height
    fig_height_px = map_height_px + mar_bottom_px + mar_top_px # [px]
   
    # Relative bottom spacing and map height
    ax_bot        = float(mar_bottom_px) / float(fig_height_px) # [-]
    ax_height     = float(map_height_px) / float(fig_height_px) # [-]
    
    # Figure size in inches for matplotlib
    fig_height_in = float(fig_height_px) / float(fig_dpi)
    fig_width_in  = float(fig_width_px ) / float(fig_dpi)

    return (ax_left, ax_bot, ax_width, ax_height, fig_height_in, fig_width_in, fig_height_px, fig_width_px)

## Function to create maps
def create_maps(olga, wrfout, dom, times):
    filter = True    # Apply Gaussian filter to variables
    fsigma = 0.5     # std dev of Gaussian filter 

    # Read WRF output
    wrf = readwrf_all(olga, wrfout, times[0], times[-1]) 

    # Get the domain properties, either from settings, or make estimation
    map_lat, map_lon, map_width, map_height = getDomain(olga, dom, wrf)

    # Calculate figure properties
    ax_left, ax_bot, ax_width, ax_height, fig_height_in, fig_width_in, fig_height_px, fig_width_px = \
        getFigureProperties(map_height, map_width)

    # Setup basemap only once and (deep)copy later for performance
    basem = Basemap(width=map_width, height=map_height,
                    rsphere=(6378137.00,6356752.3142),\
                    resolution=olga.map_res[dom], area_thresh=10., projection='lcc',\
                    lat_1=map_lat, lat_0=map_lat, lon_0=map_lon)

    # Read in the OLGA logo
    olga_logo = pl.matplotlib.image.imread(olga.olgaRoot+'include/olga_lr.png')

    # Coordinates cities # TO-DO BVS: fix
    if(False):
        if(domain==1):
            db = np.genfromtxt(olgaRoot+'include/cities_eur_d1.txt',dtype='str',delimiter=',',skip_header=0)
        elif(domain==2):
            db = np.genfromtxt(olgaRoot+'include/cities_eur_d2.txt',dtype='str',delimiter=',',skip_header=0)
        city   = db[:,0]
        citys  = db[:,1]  
        citlon = np.array(db[:,2],dtype=np.float32)  
        citlat = np.array(db[:,3],dtype=np.float32)  
        citlon,citlat = basem(citlon,citlat)

    # Loop over differnt time steps 
    for t in times:
        for var in olga.map_vars[dom]:
            fig        = pl.figure(figsize=(fig_width_in, fig_height_in), dpi=fig_dpi)
            m          = deepcopy(basem)
            axloc      = [ax_left, ax_bot, ax_width, ax_height]
            ax         = fig.add_axes(axloc)
            lon,lat    = m(wrf.lon, wrf.lat)
            minlon     = lon.min()
            maxlon     = lon.max()
            minlat     = lat.min()
            maxlat     = lat.max()
            doplot     = False

            # -------------------------------------------------
            # GENERAL  METEOROLOGY
            # -------------------------------------------------
            # -------------------------------------------------
            # sea level pressure and 10m wind
            # -------------------------------------------------
            if(var == 'slpwind'):
                title  = 'Sea level pressure (hPa) + 10m wind (kts)'
                utot   = ((wrf.U10[t,:,:]*m2k)**2. + (wrf.V10[t,:,:]*m2k)**2.)**0.5
                umax   = np.ceil((((wrf.U10[:,:,:]*m2k)**2. + (wrf.V10[:,:,:]*m2k)**2.)**0.5).max()*10.)/10.
                intv   = 4 
                levs1  = np.arange(800, 1200, 2.5)     # SLP levels
                levs2  = np.arange(0, umax, umax/50.)  # Wind levels
                pfilt  = gausf(wrf.slps[t,:,:]/100., 3.0, mode='reflect')
                cf     = False 
                barbs  = m.barbs(lon[0::intv,0::intv], lat[0::intv,0::intv],\
                                 wrf.U10[t,0::intv,0::intv]*m2k, wrf.V10[t,0::intv,0::intv]*m2k,\
                                 length=5, barbcolor='r', flagcolor='r', linewidth=0.5)
                c      = m.contour(lon, lat, pfilt, levs1, linewidths=1.5, colors='k')
                cl     = pl.clabel(c, levs1, fmt='%.1f')   
                doplot = True 

            # -------------------------------------------------
            # stream lines of wind
            # -------------------------------------------------
            if(var[:4] == 'wind'):
                if(var=='wind10m'):
                    title  = '10m wind (kts)'
                    ufilt  = gausf(wrf.U10[t,:,:], fsigma,mode='reflect')
                    vfilt  = gausf(wrf.V10[t,:,:], fsigma,mode='reflect')
                    scalefac = 25.
                elif(var=='wind1000m'):
                    title  = '1000m wind'
                    ufilt  = gausf(wrf.U1000[t,:,:], fsigma,mode='reflect')
                    vfilt  = gausf(wrf.V1000[t,:,:], fsigma,mode='reflect')
                    scalefac = 35.

                utot   = ((ufilt[:,:]*m2k)**2. + (vfilt[:,:]*m2k)**2.)**0.5
                umax   = 40. 
                levs2  = np.arange(0, umax+0.01, umax/40.)  # Wind levels
                cf     = m.contourf(lon, lat, utot, levs2, extend='both', cmap=wnd)
                units  = 'kts'
                stream = m.streamplot(lon[0,:], lat[:,0], ufilt[:,:], vfilt[:,:], density=3, linewidth=utot/scalefac, color='k')
                doplot = True 
 
            # -------------------------------------------------
            # rain
            # -------------------------------------------------
            if(var == 'rr'):
                title = 'Precipitation over last hour (mm, o=convective)'
                if(t==0):
                    rr  = wrf.rr_mp[t,:,:] + wrf.rr_con[t,:,:]
                    rrc = wrf.rr_con[t,:,:]
                else:
                    rr  = (wrf.rr_mp[t,:,:] + wrf.rr_con[t,:,:]) - (wrf.rr_mp[t-1,:,:] + wrf.rr_con[t-1,:,:])  
                    rrc = wrf.rr_con[t,:,:] - wrf.rr_con[t-1,:,:]  
                levs  = np.arange(0.1, 30.0001, 0.5)
                cf    = m.contourf(lon, lat, rr, levs, alpha=1., extend='both', cmap=rain3)
                units = 'mm'
                lim   = 0.1 ; intv  = 3
                for i in range(0, wrf.nlat, intv):
                    for j in range(0, wrf.nlon, intv):
                        if(rrc[i,j] > lim and lat[i,j] > m.ymin and lat[i,j] < m.ymax and lon[i,j] > m.xmin and lon[i,j] < m.xmax):
                            pl.text(lon[i,j], lat[i,j], 'o', size=6, ha='center', va='center', color='0.1')
                doplot = True

            # -------------------------------------------------
            # Cloud cover
            # -------------------------------------------------
            if(var == 'clouds'):
                title = 'Cloud cover fraction (o=low  -=mid  |=high)'
                levs  = np.arange(0.00, 1.001, 0.1)
                intv  = 3
                cf    = m.contourf(lon, lat, wrf.ccsum[t,:,:], levs, alpha=1., cmap=cld)
                units = 'fraction'
                lim   = 0.1


                for i in range(0, wrf.nlat, intv):
                    for j in range(0, wrf.nlon, intv):
                        if(lat[i,j] > m.ymin and lat[i,j] < m.ymax and lon[i,j] > m.xmin and lon[i,j] < m.xmax): 
                            if(wrf.cclow[t,i,j] > lim):
                                pl.text(lon[i,j], lat[i,j], 'o', size=8, ha='center', va='center', color='0.1')
                            if(wrf.ccmid[t,i,j] > lim):
                                pl.text(lon[i,j], lat[i,j], '-', size=9, ha='center', va='center', color='0.2')
                            if(wrf.cchig[t,i,j] > lim):
                                pl.text(lon[i,j], lat[i,j], '|', size=8, ha='center', va='center', color='0.1')
                doplot = True

            # -------------------------------------------------
            # Shortwave downwelling radiation
            # TO-DO: make it potential (i.e. max potential swd / actual swd)
            # -------------------------------------------------
            if(var == 'swd'):
                title = 'Shortwave incoming radiation (W/m2)'
                levs  = np.arange(0, 1000.01, 50)
                data  = gausf(wrf.swd[t,:,:], fsigma,mode='reflect') if filter else wrf.swd[t,:,:]
                cf    = m.contourf(lon, lat, data, levs, extend='both', cmap=wup)
                units = 'W/m2'
                doplot = True 
 
            # -------------------------------------------------
            # GLIDING SPECIFIC
            # -------------------------------------------------
            # -------------------------------------------------
            # Convection maps 
            # -------------------------------------------------
            if(var == 'convection'):
                title = 'asdf'
                cf = False

                cmap = ListedColormap([(205./256., 45./256., 45./256.),\
                                       (73. /256., 68./256., 68./256.),\
                                       (186./256.,186./256.,186./256.),\
                                       (119./256.,221./256.,255./256.),\
                                       (119./256.,255./256.,154./256.)], 'indexed')


                if(t==0):
                    rr  = wrf.rr_mp[t,:,:] + wrf.rr_con[t,:,:]
                    rrc = wrf.rr_con[t,:,:]
                else:
                    rr  = (wrf.rr_mp[t,:,:] + wrf.rr_con[t,:,:]) - (wrf.rr_mp[t-1,:,:] + wrf.rr_con[t-1,:,:])  
                    rrc = wrf.rr_con[t,:,:] - wrf.rr_con[t-1,:,:]  

                # 1. filter fields
                rain  = gausf(rr[:,:], fsigma, mode='reflect')
                cumul = gausf(wrf.zct[t,:,:] - wrf.zi[t,:,:], fsigma, mode='reflect')

                # Create colored background
                bg = np.zeros((wrf.nlat, wrf.nlon))
                for i in range(wrf.nlat):
                    for j in range(wrf.nlon):
                        if(rain[i,j] > 0.1): # Rain
                            bg[i,j] = 0.5
                        elif(wrf.swd_frac[t,i,j] < 0.1): 
                            bg[i,j] = 1.5
                        elif(wrf.swd_frac[t,i,j] < 0.25):
                            bg[i,j] = 2.5
                        elif(wrf.hfx[t,i,j] > 0. and cumul[i,j] < 10):
                            bg[i,j] = 3.5                  
                        elif(wrf.hfx[t,i,j] > 0. and cumul[i,j] > 10):
                            bg[i,j] = 4.5                  
                bg = np.ma.masked_where(bg==0, bg)

                cf    = m.pcolormesh(lon, lat, bg, vmin=0, vmax=5, cmap=cmap)
                units = '-'

 
                intv = 3
                for i in range(0, wrf.nlat, intv):
                    for j in range(0, wrf.nlon, intv):
                        if(lat[i,j] > m.ymin and lat[i,j] < m.ymax and lon[i,j] > m.xmin and lon[i,j] < m.xmax): 
                            if(wrf.hfx[t,i,j] > 0):
                                pl.text(lon[i,j], lat[i,j], int((wrf.zi[t,i,j]+wrf.hgt[i,j])/100), size=6, ha='center', va='center', color='0.3')
                doplot = True

            # -------------------------------------------------
            # Updraft velocity wstar
            # -------------------------------------------------
            if(var == 'wstar'):
                title = 'Updraft velocity (m/s)'
                levs  = np.arange(0, 3.001, 0.5)
                data  = gausf(wrf.wstar[t,:,:], fsigma,mode='reflect') if filter else wrf.wstar[t,:,:]
                cf    = m.contourf(lon, lat, data, levs, extend='both', cmap=wup)
                units = 'm/s'
                doplot = True 

            # -------------------------------------------------
            # TEST: z/L as indicator for cloud streets / roll convection
            # -------------------------------------------------
            if(var == 'zol'):
                title = 'z/L [-]'
                levs  = np.arange(-100., 0.01, 5.)
                data  = gausf((wrf.zi[t,:,:] / wrf.L[t,:,:]),fsigma,mode='reflect') if filter else wrf.zi[t,:,:] / wrf.L[t,:,:]
                cf    = m.contourf(lon, lat, data, levs, extend='both', cmap=cent)
                doplot = True 
 
            # -------------------------------------------------
            # Top of updraft (dry convection)
            # -------------------------------------------------
            if(var == 'zidry'):
                title = 'Updraft height [m AMSL]'
                levs  = np.arange(0, 2500.01, 300.)
                data  = gausf(wrf.zi[t,:,:] + wrf.hgt[:,:], fsigma,mode='reflect') if filter else wrf.zi[t,:,:] + wrf.hgt[:,:]
                cf    = m.contourf(lon, lat, data, levs, extend='both', cmap=wup)
                units = 'm AMSL'
                doplot = True 
    
            # -------------------------------------------------
            # Cumulus depth
            # -------------------------------------------------
            if(var == 'cudepth'):
                title = 'Cumulus depth [m]'
                levs  = np.arange(0, 4000, 400)
                data  = gausf(wrf.zct[t,:,:] - wrf.zi[t,:,:], fsigma,mode='reflect') if filter else wrf.zct[t,:,:] - wrf.zi[t,:,:]
                cf =  m.contourf(lon, lat, data, levs, extend='both', cmap=wupnl)
                units = 'm'
                doplot = True 

            # -------------------------------------------------
            # Potential flight distance per day
            # -------------------------------------------------
            if((var == 'pfd') and t < np.size(wrf.date_PFD)):
                title = 'Potential flight distance [km]'
                levs  = np.arange(0, 900.1, 100)
                data  = gausf(wrf.PFD[t,:,:], fsigma,mode='reflect') if filter else wrf.PFD[t,:,:]
                cf    =  m.contourf(lon, lat, data, levs, extend='both', cmap=wup)
                units = 'km'
                doplot = True 

            # -------------------------------------------------
            # Finish plot!
            # -------------------------------------------------
            if(doplot):
                # Plot terrain height contours
                if(False):
                    levs = arange(50, 4000.01, 50.)
                    contour(lon, lat, wrf.hgt, levs, cmap=cm.PuBu, alpha=0.5)

                # Plot cities
                if(False):
                    for i in range(len(citys)):
                        m.scatter(citlon[i],citlat[i],s=4,alpha=1.0)
                        pl.text(citlon[i],citlat[i],' '+citys[i],size=8,ha='left',va='center',color='0.1')

                # Plot airspace
                if(False):
                    plotairspace(m)

                # If filled contour, add colorbar
                if(cf != False):
                    cax = pl.axes([ax_left+ax_width+0.02, ax_bot+ 0.1*ax_height, 0.02, 0.8*ax_height])
                    cb  = pl.colorbar(cf, drawedges=False, cax=cax)
                    cb.ax.tick_params(labelsize = 8) 
                    cb.ax.yaxis.set_tick_params(width = 0)
                    cb.outline.set_color('white')
                    cb.set_label(units)
                pl.axes(ax)
          
                # Draw coast and country outlines 
                m.drawcoastlines(linewidth=1, color='0.3')
                m.drawcountries(linewidth=1, color='0.3')
                m.drawmapboundary()
            
                # Draw rivers 
                if(olga.drawRivers[dom] == True):
                    m.drawrivers(linewidth=0.5, color='#0066FF')
 
                # Add description variables, units, OLGA label, logo, etc. 
                pl.figtext(ax_left, ax_bot/2,'OLGA: %s GFS-initiated WRF-ARW forecast [olga.vanstratum.com]'\
                           %olga.map_desc[dom],size=6,ha='left', va='center')
                if(var=='pfd'):
                    subtitle = 'cumulative distance over: %s'%str(wrf.date_PFD[t]) 
                else:
                    subtitle = 'valid: %s UTC'%str(wrf.datetime[t])

                # Add figure info (variable, units, time)
                ax.set_title(title, loc='left', size=8)
                ax.set_title(subtitle, loc='right', size=8)

                # Add OLGA logo
                pl.figimage(olga_logo, fig_width_px-45, 6)

                # Get name and save figure 
                if(var=='pfd'):
                    tmp    = '%06i'%(t*24.) 
                else:
                    xtime  = wrf.time[t] / 3600.
                    hour   = int(np.floor(xtime))
                    minute = int((xtime - hour) * 60.)
                    tmp    = str(hour).zfill(4) + str(minute).zfill(2)

                #name = '%s%04i%02i%02i_t%02iz/%s_d%i_%s.png'%(olga.figRoot,olga.year,olga.month,olga.day,olga.cycle,tmp,dom+1,var)
                name = '%s%04i%02i%02i_t%02iz/%s_%02i_%s.png'%(olga.figRoot, olga.year, olga.month, olga.day, olga.cycle, var, dom+1, tmp)

                pl.savefig(name, dpi=fig_dpi)
        
            # Cleanup!
            fig.clf()
            pl.close()
            gc.collect()

