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

from readwrf import *
from colormaps import *
from tools import *
from skewtlogp import *

## Function to create maps
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

              fig=skewtlogp(sset)
          
              nameo = 'figures/'+ date + '/d' + str(domain) + '_sound' + str(stype) +'_' + name + '_' + str(t).zfill(2) + '.png'
              pl.savefig(nameo)

              # Cleanup!
              fig.clf()
              pl.close()
              gc.collect()

