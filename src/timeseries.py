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
from copy import deepcopy

from readwrf import *
from colormaps import *
from tools import *

## Function to create time series
def create_timeseries(olga,wrfout,dom,times):

    # Create 4-color colormap
    wupd = cmap_discrete(wup,np.linspace(0,1,4))
    olga_logo = pl.matplotlib.image.imread(olga.olgaRoot+'include/olga_left.png')

    # Loop over requested locations
    for name, longName, lon, lat, type in zip(olga.soundLoc[dom].shortName, olga.soundLoc[dom].longName, olga.soundLoc[dom].lon, olga.soundLoc[dom].lat, olga.soundLoc[dom].type):
        if(type == 0 or type == 2):
            # Read WRF data
	    d = readwrf_loc(olga,wrfout,lon,lat,times[0],times[-1])

            # Loop over diffent day (if available):
            for t in range(np.size(d.t0_ana)):
                t0 = d.t0_ana[t]
                t1 = d.t1_ana[t]
                x0 = d.hour[t0]
                x1 = d.hour[t1]

                fig = pl.figure(figsize=(olga.fig_width_px/float(olga.fig_dpi), olga.fig_width_px/float(olga.fig_dpi)))
                #                   L    B    R    T    ws  hs
                fig.subplots_adjust(0.10,0.11,0.98,0.88,0.35,0.55)
                pl.figtext(0.5,0.95,'%s [%.2fN, %.2fE]'%(longName,lat,lon),size=9,ha='center')
                pl.figtext(0.5,0.93,'%s'%(d.date[t0]),size=8,ha='center')
                gs = pl.matplotlib.gridspec.GridSpec(4,2,height_ratios=[1.,1.,0.4,1],width_ratios=[2,1])

                # -------------------------------------------------
                # Updraft velocity / height
                # -------------------------------------------------
                ax = pl.subplot(gs[:2,0]); modplot(ax)
                ax.set_title('Updraft velocity and height',loc='left')
                zs = d.z[0,0]  # terrain height (lowest half level)
                wm = 3.5       # scaling for colormap
                bw1 = 0.6      # width of sub-cloud updrafts
                bw2 = 0.8      # width of cloud updrafts
                for i in range(t0,t1+1):
                    if(d.wglider[i] > 0.0):
                        color = wup((np.floor(d.wglider[i])+0.5)/wm)
                        pl.bar(d.hour[i]-0.5*bw1, d.zi[i],         width=bw1, bottom=zs,         color=color,           edgecolor='none')    
                        pl.bar(d.hour[i]-0.5*bw2, d.ct[i]-d.zi[i], width=bw2, bottom=d.zi[i]+zs, color='k',  alpha=0.3, edgecolor='none')    

                # Add sort-of colorbar
                wups = ([0.5,1.5,2.5,3.5])
                names = (['0-1 m/s','1-2 m/s','2-3 m/s','>3 m/s'])
                for wu,nam in zip(wups,names):
                    pl.scatter([-10],[300],color=wup(wu/wm),label=nam)
                pl.legend(frameon=False,loc=2)  
                pl.xlim(d.hour[t0],d.hour[t1])
                pl.ylim(0,3000)
                pl.ylabel('z [m AMSL]')
                pl.xticks(np.arange(d.hour[t0],d.hour[t1]+0.001,2))
                pl.yticks(np.arange(0,3000.001,500))
                pl.grid(linestyle=':',color='0.5')

                #ax2 = ax.twinx()
                #ax2.plot(d.hour[t0:t1+1],d.cPFD[t0:t1+1],linewidth=2,color='0.5',label='CFD')
                #pl.xlim(d.hour[t0],d.hour[t1])
                #pl.ylim(0,900)
                #pl.xticks(np.arange(d.hour[t0],d.hour[t1]+0.001,2))
                #pl.legend(frameon=False,loc=1)
                #pl.ylabel('Cumulative flight distance (CFD) [km]')
                #pl.yticks(np.arange(0,900.001,150))

                # -------------------------------------------------
                # Pressure
                # -------------------------------------------------
                ax = pl.subplot(gs[0,1]); modplot(ax)
                ax.set_title('Surface pressure',loc='left')
                slp = d.slps[t0:t1]/100.
                pl.plot(d.hour[t0:t1], slp, 'k-')
                pl.xlim(d.hour[t0],d.hour[t1])
                pl.ylabel('hPa')
                pl.xticks(np.arange(d.hour[t0],d.hour[t1]+0.001,2))
                if(slp.max() - slp.min() < 10):
                    pl.ylim(roundNumber(slp.mean(), 1, DOWN)-5, roundNumber(slp.mean(), 1, UP)+5)

                # -------------------------------------------------
                # Air temperature / dew poiunt temperature
                # -------------------------------------------------
                ax = pl.subplot(gs[1,1]); modplot(ax)
                ax.set_title('T and Td at 2m',loc='left')
                pl.plot(d.hour[t0:t1], d.T2[t0:t1]-273.15, 'r-', label='T')
                pl.plot(d.hour[t0:t1], d.Td2[t0:t1]-273.15, 'b-', label='Td')
                pl.ylabel('celcius')
                pl.xlim(d.hour[t0],d.hour[t1])
                pl.xticks(np.arange(d.hour[t0],d.hour[t1]+0.001,2))
                pl.legend(frameon=False,loc=2)

                # -------------------------------------------------
                # Rain
                # -------------------------------------------------
                ax = pl.subplot(gs[2,1]); modplot(ax)
                ax.set_title('Rain',loc='left')
                pl.bar(d.hour[t0:t1], d.drr_mp[t0:t1], color='b', edgecolor='none')
                pl.bar(d.hour[t0:t1], d.drr_con[t0:t1], bottom=d.drr_mp[t0:t1], color='g', edgecolor='none')
                pl.ylabel('mm')
                pl.xlim(d.hour[t0],d.hour[t1])
                pl.xticks(np.arange(d.hour[t0],d.hour[t1]+0.001,2))
                pl.ylim(0,max(1,roundNumber((d.drr_mp[t0:t1]+d.drr_con[t0:t1]).max(), 1, UP))) 

                # -------------------------------------------------
                # Shortwave radiation
                # -------------------------------------------------
                ax = pl.subplot(gs[3,1]); modplot(ax)
                ax.set_title('Shortwave radiation',loc='left')
                pl.plot(d.hour[t0:t1], d.swd[t0:t1], 'k-')
                pl.plot(d.hour[t0:t1], d.swdc[t0:t1], 'k:')
                pl.ylabel('W/m2')
                pl.xlabel('time UTC [h]')
                pl.xlim(d.hour[t0],d.hour[t1])
                pl.xticks(np.arange(d.hour[t0],d.hour[t1]+0.001,2))

                # -------------------------------------------------
                # Cloud cover
                # -------------------------------------------------
                ax = pl.subplot(gs[2,0]); modplot(ax)
                ax.set_title('Cloud cover',loc='left')
                pl.pcolormesh(d.ccl,cmap=pl.cm.bone_r,vmin=0,vmax=1)
                pl.text(d.hour[t0]-0.4,0.5,'low',size=7,ha='right',va='center')
                pl.text(d.hour[t0]-0.4,1.5,'middle',size=7,ha='right',va='center')
                pl.text(d.hour[t0]-0.4,2.5,'high',size=7,ha='right',va='center')
                pl.xlim(d.hour[t0],d.hour[t1])
                ax.set_yticks([])
                pl.xticks(np.arange(d.hour[t0],d.hour[t1]+0.001,2))

                #ax2 = ax.twinx()
                #ax2.plot(d.hour[t0:t1+1],d.swd_frac[t0:t1+1]*100.,linewidth=2,color='y',label='% sun')
                #pl.xlim(d.hour[t0],d.hour[t1])
                #pl.ylim(0,105)
                #pl.yticks([0,25,50,75,100])
                #pl.legend(frameon=False,loc=4)
                #pl.ylabel('Sun [%]')

                # -------------------------------------------------
                # Wind
                # -------------------------------------------------
                ax = pl.subplot(gs[3,0]); modplot(ax)
                ax.set_title('Wind',loc='left')
                k500  = key_nearest(d.zf[0,:],500)
                k1000 = key_nearest(d.zf[0,:],1000)
                k2000 = key_nearest(d.zf[0,:],2000)
                pl.barbs(d.hour[t0:t1+1],0.5,d.u10[t0:t1+1]      *m2k,d.v10[t0:t1+1]      *m2k,length=5,linewidth=0.5,pivot='middle')  
                pl.barbs(d.hour[t0:t1+1],1.5,d.u  [t0:t1+1,k500] *m2k,d.v  [t0:t1+1,k500] *m2k,length=5,linewidth=0.5,pivot='middle')  
                pl.barbs(d.hour[t0:t1+1],2.5,d.u  [t0:t1+1,k1000]*m2k,d.v  [t0:t1+1,k1000]*m2k,length=5,linewidth=0.5,pivot='middle')  
                pl.barbs(d.hour[t0:t1+1],3.5,d.u  [t0:t1+1,k2000]*m2k,d.v  [t0:t1+1,k2000]*m2k,length=5,linewidth=0.5,pivot='middle')  
                pl.text(d.hour[t0]-0.4,0.5,'10m',  size=7,ha='right',va='center')
                pl.text(d.hour[t0]-0.4,1.5,'500m', size=7,ha='right',va='center')
                pl.text(d.hour[t0]-0.4,2.5,'1000m',size=7,ha='right',va='center')
                pl.text(d.hour[t0]-0.4,3.5,'2000m',size=7,ha='right',va='center')
                pl.xlim(0,24)
                pl.xlim(d.hour[t0],d.hour[t1])
                ax.set_yticks([])
                pl.xticks(np.arange(d.hour[t0],d.hour[t1]+0.001,2))
                pl.xlabel('time UTC [h]')

                # Add logo (105px wide) 
                pl.figimage(olga_logo, 10, olga.fig_width_px-45)
                pl.figtext(0.01, 0.011, '%s'%(olga.map_desc[dom]), size=7, ha='left')

                # Save figure
                tmp  = '%06i'%(olga.islice*24.)
                #name = '%s%04i%02i%02i_t%02iz/%s_d%i_tser_%s.png'%(olga.figRoot,olga.year,olga.month,olga.day,olga.cycle,tmp,dom+1,name)
                name = '%s%04i%02i%02i_t%02iz/tser_%s_%02i_%s.png'%(olga.figRoot, olga.year, olga.month, olga.day, olga.cycle, name, dom+1, tmp)
                pl.savefig(name, dpi=olga.fig_dpi)
