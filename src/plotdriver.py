#
# Copyright (c) 2013-2016 Bart van Stratum
# Copyright (c) 2015-2016 Roel Baardman
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

# Use non-graphical backend to prevent figures popping up:
import matplotlib as mpl
mpl.use('Agg') 

import sys
import os
import glob
import numpy as np
import pylab as pl

pl.matplotlib.rc('font', size=9)
pl.matplotlib.rc('legend', fontsize=8)

from maps import *         # plot maps
from timeseries import *   # plot time series
from sounding import *     # plot soundings

def plot_driver(olga,dom,ptype):
    # Check if wrfout file is available
    wrfout = '%s%04i%02i%02i_t%02iz_d%i.nc'%(olga.wrfDataRoot,olga.year,olga.month,olga.day,olga.cycle,dom+1)
    if(not os.path.exists(wrfout)):
        sys.exit('Cant find %s'%wrfout)

    # Create image directory if not available
    figpath = '%s%04i%02i%02i_t%02iz'%(olga.figRoot,olga.year,olga.month,olga.day,olga.cycle)
    if not os.path.exists(figpath):
        try:
            os.makedirs(figpath)
        except:
            print 'BvS have to solve this for parallel execution...' # Jul2014: should be solved now..

    # Create array with times to plot:
    t0=olga.t0/(olga.dt_output[dom]/60.)
    t1=olga.t1/(olga.dt_output[dom]/60.)
    times = np.arange(t0,t1+1e-9,dtype=np.int)

    if(ptype == "maps"):
        createMaps(olga,wrfout,dom,times) 
    elif(ptype=="time"):
        create_timeseries(olga,wrfout,dom,times)
    elif(ptype=="sounding"):
        create_sounding(olga,wrfout,dom,times)
    else:
        sys.exit('type=%s not supported'%type)

