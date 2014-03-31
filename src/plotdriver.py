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

# Use non-graphical backend to prevent figures popping up:
import matplotlib as mpl
mpl.use('Agg') 

import sys
import os
import glob
import numpy as np
import pylab as pl
from netCDF4 import Dataset 

from settings import *     # OLGA settings
from maps import *         # plot maps
from timeseries import *   # plot time series
from sounding import *     # plot soundings


"""
Check if requested WRF file is available
"""
def checkWRFoutput(year,month,day,nrun,cycle,dom,date):
    wrfout = wrfDataRoot+date+'_t'+str(cycle).zfill(2)+'z_d'+str(dom)+'_'+str(nrun)+'.nc'

    if(not os.path.exists(wrfout)):
        sys.exit('cant find %s'%wrfout)
    else:
        nt = np.size(Dataset(wrfout,'r').variables["XTIME"][:])
        return (nt,wrfout)

"""
Create maps
"""
def maps(wrfout,dom,date,times,nrun):
    # Variables to plot per domain:
    if(dom==1):
        variables = (['clouds','rr','wind10m','wind1000m']) 
    elif(dom==2):
        variables = (['pfd','wstar','zidry','cudepth']) 
    create_maps(wrfout,dom,date,times,variables,nrun,filter=True)

"""
Create meteograms
"""
def timeseries(wrfout,dom,date,times,nruns):
    # Get list of locations
    l = np.genfromtxt(olgaRoot+'include/analysis_locs.txt', delimiter=',', dtype=("|S10",float,float))
    names = l["f0"]
    lons =  l["f1"]
    lats =  l["f2"]
    create_tser(wrfout,dom,date,names,lons,lats,nruns)

"""
Create soundings
"""
def soundings(wrfout,dom,date,times):
    # Get list of locations
    l = np.genfromtxt(olgaRoot+'include/sounding_locs.txt', delimiter=',', dtype=("|S15",float,float))
    names = l["f0"]
    lons =  l["f1"]
    lats =  l["f2"]
    create_sounding(wrfout,dom,date,times,names,lons,lats)

def makeplots(year,month,day,nrun,cycle,dom,type):

    date = str(year)+str(month).zfill(2)+str(day).zfill(2)

    # Check if wrfout available, and get file path and number of output time steps
    nt,wrfout = checkWRFoutput(year,month,day,nrun,cycle,dom,date)

    # Create image directory if not available
    if not os.path.exists(figRoot+'/'+date):
        try:
            os.makedirs(figRoot+'/'+date)
        except:
            print 'BvS have to solve this for parallel execution...'

    # Which times to plot?
    times = np.arange(0,nt,1)

    if(type == "maps"):
        maps(wrfout,dom,date,times,nrun) 
    elif(type=="time"):
        timeseries(wrfout,dom,date,times,nrun) 
    elif(type=="sounding"):
        soundings(wrfout,dom,date,times) 
    else:
        sys.exit('type=%s not supported'%type)

"""
plotdriver can either be called directly
in which case this code is triggered
or called from e.g. WRFdriver with doplots()
"""
if __name__ == "__main__":
    # BvS move to matplotlibrc
    pl.matplotlib.rc('font', size=9)
    pl.matplotlib.rc('legend', fontsize=8)

    # Get command line arguments
    if(len(sys.argv) != 8):
        sys.exit('provide input: "YYYY MM DD run cycle dom type}"')
    else:
        year   = int(sys.argv[1])
        month  = int(sys.argv[2])
        day    = int(sys.argv[3])
        nrun   = int(sys.argv[4])
        cycle  = int(sys.argv[5])
        dom    = int(sys.argv[6])
        type   = sys.argv[7]

    makeplots(year,month,day,nrun,cycle,dom,type)


