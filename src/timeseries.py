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
    olga_logo = pl.matplotlib.image.imread(olga.olgaRoot+'include/olga_lr.png')

    # Loop over requested locations
    for name,lon,lat in zip(olga.meteogr_name[dom][:],olga.meteogr_lon[dom][:],olga.meteogr_lat[dom][:]):
        # Read WRF data
	d = readwrf_loc(olga,wrfout,lon,lat,times[0],times[-1])

        # Loop over diffent day (if available):
        for t in range(np.size(d.t0_ana)):
            t0 = d.t0_ana[t]
            t1 = d.t1_ana[t]

            fig = pl.figure(figsize=(6.5,6.0))
            #                   L    B    R    T    ws  hs
            fig.subplots_adjust(0.10,0.11,0.90,0.88,0.36,0.36)
            pl.figtext(0.5,0.95,'%s [%.2fN, %.2fE]'%(name,lat,lon),size=9,ha='center')
            pl.figtext(0.5,0.93,'%s'%(d.date[t0]),size=8,ha='center')
            gs = pl.matplotlib.gridspec.GridSpec(3,1,height_ratios=[4,1,1.5])

            # -------------------------------------------------
            # Updraft velocity / height
            # -------------------------------------------------
            ax = pl.subplot(gs[0])
            ax.set_title('Updraft velocity and height',loc='left')
            zs = d.z[0,0]  # terrain height (lowest half level)
            wm = 3.5       # scaling for colormap
            for i in range(t0,t1+1):
                if(d.wstar[i] > 0.0):
                    pl.bar(d.hour[i]-0.2,d.zi[i],width=0.4,bottom=zs,color=wup((np.floor(d.wstar[i])+0.5)/wm),edgecolor='none')    
                    pl.bar(d.hour[i]-0.225,d.ct[i]-d.zi[i],width=0.45,bottom=d.zi[i]+zs,color='k',alpha=0.3,edgecolor='none')    
            # Add sort-of colorbar
            wups = ([0.5,1.5,2.5,3.5])
            names = (['0-1 m/s','1-2 m/s','2-3 m/s','>3 m/s'])
            for wu,nam in zip(wups,names):
                pl.scatter([-10],[300],color=wup(wu/wm),label=nam)
            pl.legend(frameon=False,loc=2)  
            pl.xlim(d.hour[t0],d.hour[t1])
            pl.ylim(0,3000)
            modplot(ax,removeaxis=['top']) 
            pl.ylabel('z [m AMSL]')
            pl.xticks(np.arange(d.hour[t0],d.hour[t1]+0.001,1))

            ax2 = ax.twinx()
            ax2.plot(d.hour[t0:t1+1],d.cPFD[t0:t1+1],linewidth=2,color='0.5',label='CFD')
            pl.xlim(d.hour[t0],d.hour[t1])
            pl.ylim(0,1000)
            pl.legend(frameon=False,loc=1)
            pl.ylabel('Cumulative flight distance (CFD) [km]')

            # -------------------------------------------------
            # Cloud cover
            # -------------------------------------------------
            ax = pl.subplot(gs[1])
            ax.set_title('Cloud cover',loc='left')
            pl.pcolormesh(d.ccl,cmap=pl.cm.bone_r,vmin=0,vmax=1)
            pl.text(d.hour[t0]-0.4,0.5,'low',size=7,ha='right',va='center')
            pl.text(d.hour[t0]-0.4,1.5,'middle',size=7,ha='right',va='center')
            pl.text(d.hour[t0]-0.4,2.5,'high',size=7,ha='right',va='center')
            pl.xlim(d.hour[t0],d.hour[t1])
            pl.ylabel('z [m]')
            modplot(ax,removeaxis=['top']) 
            ax.set_yticks([])
            pl.xticks(np.arange(d.hour[t0],d.hour[t1]+0.001,1))

            ax2 = ax.twinx()
            ax2.plot(d.hour[t0:t1+1],d.swd_frac[t0:t1+1]*100.,linewidth=2,color='y',label='% sun')
            pl.xlim(d.hour[t0],d.hour[t1])
            pl.ylim(0,105)
            pl.yticks([0,25,50,75,100])
            pl.legend(frameon=False,loc=4)
            pl.ylabel('Sun [%]')

            # -------------------------------------------------
            # Wind
            # -------------------------------------------------
            ax = pl.subplot(gs[2])
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
            modplot(ax)
            ax.set_yticks([])
            pl.xticks(np.arange(d.hour[t0],d.hour[t1]+0.001,1))
            pl.xlabel('time UTC [h]')

            # Add logo 
            w=650;h=600
            pl.figimage(olga_logo,10,6)
            pl.figtext(0.5,0.010,'OLGA: 6 x 6 km GFS-initiated WRF-ARW forecast [olga.vanstratum.com]',size=7,ha='center')

            # Save figure
            tmp  = '%04i'%(olga.islice*24.)
            name = '%s%04i%02i%02i_t%02iz/%s_d%i_tser_%s.png'%(olga.figRoot,olga.year,olga.month,olga.day,olga.cycle,tmp,dom+1,name)
            pl.savefig(name)

