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
import numpy as np
import pylab as pl

from maps import *         # plot maps
from timeseries import *   # plot time series
from sounding import *     # plot soundings

## System settings
wrfdataroot = '../dataWRF' # path to WRF output
imagedir    = 'figures'    # path to put figures

# BvS move to matplotlibrc
pl.matplotlib.rc('font', size=9)
pl.matplotlib.rc('legend', fontsize=8)

## Get command line arguments
if(len(sys.argv) != 6):
    sys.exit('provide input: "YYYY MM DD {1/2} {maps/tser/sounding}" (year month day domain type)')
else:
    year   = int(sys.argv[1])
    month  = int(sys.argv[2])
    day    = int(sys.argv[3])
    dom    = int(sys.argv[4])
    type   = sys.argv[5]
    dt     = 1

## Check if wrfout available
date   = str(year)+str(month).zfill(2)+str(day).zfill(2)
wrfout = wrfdataroot+'/'+date+'/wrfout_d0'+str(dom)+'_'+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_00:00:00' 
if(not os.path.isfile(wrfout)):
    sys.exit('file %s not available'%wrfout) 

## Create image directory if not available
if not os.path.exists(imagedir+'/'+date):
    os.makedirs(imagedir+'/'+date)

## Create maps
if(type == "maps"):
    # Variables to plot per domain:
    if(dom==1):
        variables = (['clouds','rr','slpwind','clouds2']) 
    elif(dom==2):
        variables = (['pfd','wstar','zidry','cudepth']) 
    create_maps(wrfout,dom,date,0,24,dt,variables,filter=True)

## Create meteograms
elif(type=="tser"):
    # Get list of locations
    l = np.genfromtxt('analysis_locs.txt', delimiter=',', dtype=("|S10",float,float))
    names = l["f0"]
    lons =  l["f1"]
    lats =  l["f2"]
    create_tser(wrfout,dom,date,names,lons,lats)

## Create soundings
elif(type=="sounding"):
    # Get list of locations
    l = np.genfromtxt('sounding_locs.txt', delimiter=',', dtype=("|S15",float,float))
    names = l["f0"]
    lons =  l["f1"]
    lats =  l["f2"]
    times = np.arange(0,25,1)
    create_sounding(wrfout,dom,date,names,lons,lats,times)

## Exception
else:
    sys.exit('type=%s not supported'%type)
