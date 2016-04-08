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
from scipy.ndimage.filters import gaussian_filter as gaussianFilter
from multiprocessing import Process

# Custom modules
from readwrf import *
from colormaps import *
from tools import *

smoothPlot = True  # Apply Gaussian filter to variables
#fsigma = 0.5     # std dev of Gaussian filter 

""" 
* Set the map properties, either from
* settings, or estimate from WRF output 
""" 
def getDomain(olga, dom, wrf):
    if(olga.map_lat[dom] == -1):
        map_lat = np.average(wrf.lat)
        print('latitude domain %i = %f'%(dom, map_lat))
    else:
        map_lat = olga.map_lat[dom]

    if(olga.map_lon[dom] == -1):
        map_lon = np.average(wrf.lon)
        print('longitude domain %i = %f'%(dom, map_lon))
    else:
        map_lon = olga.map_lon[dom]

    if(olga.map_width[dom] == -1):
        map_width =  1.05 * haversine(wrf.lon.min(), map_lat, wrf.lon.max(), map_lat) * 1000
        print('width domain %i = %f'%(dom, map_width))
    else:
        map_width = olga.map_width[dom]

    if(olga.map_height[dom] == -1):
        map_height =  1.05 * haversine(map_lon, wrf.lat.min(), map_lon, wrf.lat.max()) * 1000
        print('height domain %i = %f'%(dom, map_height))
    else:
        map_height = olga.map_height[dom]

    return (map_lat, map_lon, map_width, map_height)

""" 
* Calculate the figure properties to ensure:
* 1. Prescribed figure width in pixels
* 2. Fixed margins for labels, colorbar, etc 
"""
def getFigureProperties(olga, map_height, map_width):
    # Aspect ratio of map
    map_aspect    = (float(map_height) / float(map_width))
   
    # Absolute size map
    map_width_px  = olga.fig_width_px - olga.map_left_px - olga.map_right_px # [px]
    map_height_px = map_width_px * map_aspect # [px]

    # Relative left spacing and map width
    ax_left       = float(olga.map_left_px)  / float(olga.fig_width_px) 
    ax_width      = float(map_width_px) / float(olga.fig_width_px)
   
    # Figure height
    fig_height_px = map_height_px + olga.map_bottom_px + olga.map_top_px # [px]
   
    # Relative bottom spacing and map height
    ax_bot        = float(olga.map_bottom_px) / float(fig_height_px) # [-]
    ax_height     = float(map_height_px) / float(fig_height_px) # [-]
    
    # Figure size in inches for matplotlib
    fig_height_in = float(fig_height_px) / float(olga.fig_dpi)
    fig_width_in  = float(olga.fig_width_px ) / float(olga.fig_dpi)

    return (ax_left, ax_bot, ax_width, ax_height, fig_height_in, fig_width_in, fig_height_px)

""" 
* Create and save a single figure
"""
def createFigure(olga, dom, wrf, basem, var, t, figwi, fighi, dpi, axl, axb, axw, axh):
    # This is ugly :( TODO: fix
    if(var == 'pfd' or var == 'pfd2'):
        nFig = np.size(olga.pfdNames)
    else:
        nFig = 1

    for n in range(nFig):
        fig        = pl.figure(figsize=(figwi, fighi), dpi=dpi)
        m          = deepcopy(basem)
        axloc      = [axl, axb, axw, axh]
        ax         = fig.add_axes(axloc)
        lon,lat    = m(wrf.lon, wrf.lat)
        minlon     = lon.min()
        maxlon     = lon.max()
        minlat     = lat.min()
        maxlat     = lat.max()
        doplot     = False

        # Default color line country boundaries
        countryLines = '0.3'

        # -------------------------------------------------
        # stream lines of wind
        # -------------------------------------------------
        if(var[:4] == 'wind'):
            doplot = True 
            units  = 'kts'
            fsigma = 0.5

            if(var=='wind10m'):
                title = '10m wind (kts)'
                ufilt = gaussianFilter(wrf.U10[t,:,:], fsigma, mode='reflect')
                vfilt = gaussianFilter(wrf.V10[t,:,:], fsigma, mode='reflect')
                sf    = 25.
            elif(var=='wind1000m'):
                title = '1000m wind'
                ufilt = gaussianFilter(wrf.U1000[t,:,:], fsigma, mode='reflect')
                vfilt = gaussianFilter(wrf.V1000[t,:,:], fsigma, mode='reflect')
                sf    = 35.

            umax   = 40. 
            levs   = np.arange(0, umax+0.01, 1.)
            utot   = ((ufilt[:,:]*m2k)**2. + (vfilt[:,:]*m2k)**2.)**0.5
            cf     = m.contourf(lon, lat, utot, levs, extend='both', cmap=wnd)
            stream = m.streamplot(lon[0,:], lat[:,0], ufilt[:,:], vfilt[:,:], density=3, linewidth=utot/sf, color='k')
 
        # -------------------------------------------------
        # rain
        # -------------------------------------------------
        elif(var == 'rain'):
            doplot = True
            title  = 'Precip. between t=%02i-%02i UTC (mm, //=convective)'%(wrf.hour[t-1],wrf.hour[t])
            units  = 'mm'
            fSigma = 0.5

            # At first output time instantaneous precipitation field is given
            # At consecutive steps, precipitation is accumulated. Calculate difference
            # between output time steps, so plotted field is precipitation over last period!
            if(t==0):
                rr  = wrf.rr_mp[t,:,:] + wrf.rr_con[t,:,:]
                rrc = wrf.rr_con[t,:,:]
            else:
                rr  = wrf.rr_mp[t,:,:] + wrf.rr_con[t,:,:] - wrf.rr_mp[t-1,:,:] - wrf.rr_con[t-1,:,:]  
                rrc = wrf.rr_con[t,:,:] - wrf.rr_con[t-1,:,:] 

            if(smoothPlot):
                rr  = gaussianFilter(rr,  fSigma, mode='reflect') 
                rrc = gaussianFilter(rrc, fSigma, mode='reflect') 

            # Draw shaded contours for precipitation intensity
            levs   = [0,0.1,0.2,0.4,0.8,1,2,3,4,5,6,7,8,9,10,12,14]
            cf     = m.contourf(lon, lat, rr, levs, alpha=1., extend='both', cmap=rain3nl)

            # Hatch area where precipitation is convective
            levs   = [0.1,100] 
            cf2    = m.contourf(lon, lat, rrc, levs, colors='none', hatches=['//']) 

        # -------------------------------------------------
        # Cloud cover
        # -------------------------------------------------
        elif(var == 'clouds'):
            doplot = True
            title  = 'Cloud fraction (- low  / mid  . high)'
            units  = 'fraction'
            fSigma = 0.5

            cclow  = gaussianFilter(wrf.cclow[t,:,:], fSigma, mode='reflect') if smoothPlot else wrf.cclow[t,:,:]
            ccmid  = gaussianFilter(wrf.ccmid[t,:,:], fSigma, mode='reflect') if smoothPlot else wrf.ccmid[t,:,:]
            cchig  = gaussianFilter(wrf.cchig[t,:,:], fSigma, mode='reflect') if smoothPlot else wrf.cchig[t,:,:]
            ccsum  = gaussianFilter(wrf.ccsum[t,:,:], fSigma, mode='reflect') if smoothPlot else wrf.ccsum[t,:,:]

            # Draw shaded contours for cloud cover
            levs   = np.arange(0.00, 1.001, 0.1)
            levs[0] = 0.02
            cf     = m.contourf(lon, lat, ccsum, levs, alpha=1., cmap=cm.GnBu)

            # Hatch areas of low/mid/high clouds
            levs   = [0.05,100] 
            cf2    = m.contourf(lon, lat, cclow, levs, colors='none', hatches=['-']) 
            cf2    = m.contourf(lon, lat, ccmid, levs, colors='none', hatches=['/']) 
            cf2    = m.contourf(lon, lat, cchig, levs, colors='none', hatches=['.']) 

        # -------------------------------------------------
        # Fraction of potential incoming shortwave radiation at surface
        # -------------------------------------------------
        elif(var == 'swd'):
            doplot = True
            title  = 'Incoming solar radiation (// = night)'
            units  = 'percentage'
            fSigma = 0.5

            data  = gaussianFilter(wrf.swdf[t,:,:], fSigma, mode='reflect') if smoothPlot else wrf.swdf[t,:,:]

            # Draw shaded contours for cloud cover
            levs   = np.arange(0.00, 100.01, 5)
            cf     = m.contourf(lon, lat, data*100., levs, extend='both', alpha=1., cmap=cloud)

            # Hatch night area
            levs   = [-2,-0.5] 
            cf2    = m.contourf(lon, lat, data, levs, colors='none', hatches=['/']) 

        # -------------------------------------------------
        # Updraft velocity wstar
        # -------------------------------------------------
        elif(var == 'wglider'):
            doplot = True 
            title  = 'Updraft velocity w* - %.1f m/s'%olga.sinkGlider
            units  = 'm/s'
            fSigma = 1.

            levs   = np.arange(0, 3.001, 0.5)
            data   = gaussianFilter(wrf.wglider[t,:,:], fSigma, mode='reflect') if smoothPlot else wrf.wglider[t,:,:]
            cf     = m.contourf(lon, lat, data, levs, extend='both', cmap=wup)

        # -------------------------------------------------
        # Updraft velocity TEMF
        # -------------------------------------------------
        elif(var == 'wgliderTEMF'):
            doplot = True 
            title  = 'Updraft velocity TEMF - %.1f m/s'%olga.sinkGlider
            units  = 'm/s'
            fSigma = 1.

            levs   = np.arange(0, 3.001, 0.5)
            data   = gaussianFilter(wrf.wglider2[t,:,:], fSigma, mode='reflect') if smoothPlot else wrf.wglider2[t,:,:]
            cf     = m.contourf(lon, lat, data, levs, extend='both', cmap=wup)

        # -------------------------------------------------
        # Top of updraft (dry convection)
        # -------------------------------------------------
        elif(var == 'zidry'):
            doplot = True 
            units  = 'm AMSL'
            title  = 'Updraft height [m AMSL]'
            fSigma = 1.

            levs   = np.arange(0, 2500.01, 300.)
            data   = gaussianFilter(wrf.zi[t,:,:] + wrf.hgt[:,:], fSigma, mode='reflect') if smoothPlot else wrf.zi[t,:,:] + wrf.hgt[:,:]
            cf     = m.contourf(lon, lat, data, levs, extend='both', cmap=wup)
       
        # -------------------------------------------------
        # Top of updraft (minus glider sink)
        # -------------------------------------------------
        elif(var == 'ziglider'):
            doplot = True 
            units  = 'm AMSL'
            title  = 'Updraft height with %.1f m/s sink [m AMSL]'%olga.sinkGlider
            fSigma = 1.

            levs   = np.arange(0, 2500.01, 300.)
            data   = gaussianFilter(wrf.zi2[t,:,:] + wrf.hgt[:,:], fSigma, mode='reflect') if smoothPlot else wrf.zi2[t,:,:] + wrf.hgt[:,:]
            cf     = m.contourf(lon, lat, data, levs, extend='both', cmap=wup)
 
        # -------------------------------------------------
        # Cumulus depth
        # -------------------------------------------------
        elif(var == 'cudepth'):
            doplot = True 
            title  = 'Cumulus depth [m]'
            units  = 'm'
            fSigma = 1.

            levs   = np.arange(0, 4000, 400)
            data   = gaussianFilter(wrf.zct[t,:,:] - wrf.zi[t,:,:], fSigma, mode='reflect') if filter else wrf.zct[t,:,:] - wrf.zi[t,:,:]
            cf     =  m.contourf(lon, lat, data, levs, extend='both', cmap=wupnl)

        # -------------------------------------------------
        # Potential flight distance per day
        # -------------------------------------------------
        elif((var == 'pfd') and t < np.size(wrf.date_PFD)):
            doplot = True 
            title  = 'PFD (w*) %s [km]'%olga.pfdNotes[n]
            units  = 'km'
            fSigma = 1.

            levs   = np.arange(0, 800.1, 100)
            data   = gaussianFilter(wrf.PFD1[t,n,:,:], fSigma, mode='reflect') if filter else wrf.PFD1[t,n,:,:]
            cf     =  m.contourf(lon, lat, data, levs, extend='both', cmap=wup)

        # -------------------------------------------------
        # Potential flight distance per day: TEMF corrected updraft velocity / height
        # -------------------------------------------------
        elif((var == 'pfd2') and t < np.size(wrf.date_PFD)):
            doplot = True 
            title  = 'PFD (TEMF) %s [km]'%olga.pfdNotes[n]
            units  = 'km'
            fSigma = 1.

            levs   = np.arange(0, 800.1, 100)
            data   = gaussianFilter(wrf.PFD2[t,n,:,:], fSigma, mode='reflect') if filter else wrf.PFD2[t,n,:,:]
            cf     =  m.contourf(lon, lat, data, levs, extend='both', cmap=wup)

        # -------------------------------------------------
        # Finish plot!
        # -------------------------------------------------
        if(doplot):
            # Plot terrain height contours
            if(False):
                levs = arange(50, 4000.01, 50.)
                contour(lon, lat, wrf.hgt, levs, cmap=cm.PuBu, alpha=0.5)

            # Plot cities
            if(olga.drawCities[dom]):
                for i in range(len(olga.cityLonM)):
                    m.scatter(olga.cityLonM[i], olga.cityLatM[i], s=4, alpha=1.0)
                    pl.text(olga.cityLonM[i], olga.cityLatM[i],' '+olga.cityLoc[dom].shortName[i], size=8, ha='left', va='center', color='0.1')

            # Plot airspace
            if(False):
                plotairspace(m)

            # If filled contour, add colorbar
            if(cf != False):
                cax = pl.axes([axl+axw+0.02, axb+ 0.1*axh, 0.02, 0.8*axh])
                cb  = pl.colorbar(cf, drawedges=True, cax=cax)
                cb.ax.tick_params(labelsize = 8) 
                cb.ax.yaxis.set_tick_params(width = 0)
#                cb.outline.set_color('white')
                cb.set_label(units)
            pl.axes(ax)
        
            # Draw coast and country outlines 
            m.drawcoastlines(linewidth=1.5, color=countryLines)
            m.drawcountries(linewidth=1.5, color=countryLines)
        
            # Draw rivers 
            if(olga.drawRivers[dom]):
                m.drawrivers(linewidth=0.5, color='#0066FF')
 
            # Add description variables, units, OLGA label, logo, etc. 
            pl.figtext(axl, axb/2,'%s'%olga.map_desc[dom],size=7,ha='left', va='center')

            if(var=='pfd'):
                subtitle = 'cumulative distance over: %s'%str(wrf.date_PFD[t]) 
            else:
                subtitle = 'valid: %s UTC'%str(wrf.datetime[t])

            # Add figure info (variable, units, time)
            ax.set_title(title, loc='left', size=8)
            ax.set_title(subtitle, loc='right', size=8)

            # Add OLGA logo
            pl.figimage(olga.olga_logo, olga.fig_width_px-121, 4)

            # Get name and save figure 
            if(var=='pfd' or var == 'pfd2'):
                tmp    = '%06i'%(olga.islice*24*100) 
            else:
                xtime  = wrf.time[t] / 3600.
                hour   = int(np.floor(xtime))
                minute = int((xtime - hour) * 60.)
                tmp    = str(hour).zfill(4) + str(minute).zfill(2)

            if(var == 'pfd' or var == 'pfd2'):
                name = '%s%04i%02i%02i_t%02iz/%s_%s_%02i_%s.png'%(olga.figRoot, olga.year, olga.month, olga.day, olga.cycle, var, olga.pfdNames[n], dom+1, tmp)
            else:
                name = '%s%04i%02i%02i_t%02iz/%s_%02i_%s.png'%(olga.figRoot, olga.year, olga.month, olga.day, olga.cycle, var, dom+1, tmp)

            pl.savefig(name, dpi=olga.fig_dpi)
        
        # Cleanup!
        fig.clf()
        pl.close()
        gc.collect()

""" 
* This function is called from plotdriver.py, it reads the output data from WRF,
* sets up the base map (map projection, etc.) and creates the different maps 
* in parallel, to speed up execution
"""
def createMaps(olga, wrfout, dom, times):
    # Read WRF output
    wrf = readwrf_all(olga, wrfout, times[0], times[-1]) 

    # Get the domain properties, either from settings, or make estimation
    map_lat, map_lon, map_width, map_height = getDomain(olga, dom, wrf)

    # Calculate figure properties
    ax_left, ax_bot, ax_width, ax_height, fig_height_in, fig_width_in, fig_height_px = \
        getFigureProperties(olga, map_height, map_width)

    # Setup basemap only once and (deep)copy later for performance
    basem = Basemap(width=map_width, height=map_height,
                    rsphere=(6378137.00,6356752.3142),\
                    resolution=olga.map_res[dom], area_thresh=10., projection='lcc',\
                    lat_1=map_lat, lat_0=map_lat, lon_0=map_lon)

    # Read in the OLGA logo
    olga.olga_logo = pl.matplotlib.image.imread(olga.olgaRoot+'include/olga_right.png')

    # Coordinates cities/airfields/etc
    if(olga.drawCities[dom]):
        olga.cityLonM, olga.cityLatM = basem(olga.cityLoc[dom].lon, olga.cityLoc[dom].lat)

    # Serial version (13m33s):
    #for t in times:
    #    processList = []
    #    for var in olga.map_vars[dom]:
    #        createFigure(olga, dom, wrf, basem, var, t, fig_width_in, fig_height_in, olga.fig_dpi, ax_left, ax_bot, ax_width, ax_height)

    # Parallel version (4m10s): 
    for var in olga.map_vars[dom]:
        processList = []
        # For each variable, plot all maps in parallel:
        for t in range(wrf.nt):
            processList.append(Process(target=createFigure, \
                args=(olga, dom, wrf, basem, var, t, fig_width_in, fig_height_in, olga.fig_dpi, ax_left, ax_bot, ax_width, ax_height,)))
            processList[-1].start()
        # Join all processes:
        for process in processList:
            process.join()
