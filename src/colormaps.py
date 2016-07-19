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

import numpy as np
import pylab as pl
import glob
import sys
import os
from matplotlib.colors import LinearSegmentedColormap, ColorConverter

def cmap_discrete(cmap_in,ints):
    return pl.get_cmap(cmap_in)(ints)

# A bit hacky....
filedir = os.path.dirname(__file__)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create colormap from NCL's .rgb files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def cmap_ncl(name):
    colormap = os.path.join(filedir, '../include/colormaps/%s.rgb'%name)

    if not os.path.isfile(colormap):
        sys.exit('cant find file %s.rgb! exit'%name) 
      
    cols = []
    f = open(colormap, 'r')
    l = 0
    for line in f:
        if(l==0):
            ncolors = int(line.split('=')[-1])
        elif(line.split()[0]!='#'):
            for i in range(3):
                cols.append(int(line.split()[i])) 
        l+=1

    R    = np.array(cols[0::3])/256.
    G    = np.array(cols[1::3])/256.
    B    = np.array(cols[2::3])/256.

    x0 = np.linspace(0,1,ncolors)
    cmap_dict = {}
    cmap_dict['red'] = [(x0[i],R[i],R[i]) for i in range(len(R))]
    cmap_dict['green'] = [(x0[i],G[i],G[i]) for i in range(len(G))]
    cmap_dict['blue'] = [(x0[i],B[i],B[i]) for i in range(len(B))]
    mymap = LinearSegmentedColormap('mymap',cmap_dict)
    return mymap

def make_colormap(colors):
    from matplotlib.colors import LinearSegmentedColormap, ColorConverter
    from numpy import sort

    z  = np.array(sorted(colors.keys()))
    n  = len(z)
    z1 = min(z)
    zn = max(z)
    x0 = (z - z1) / (zn - z1)

    CC = ColorConverter()
    R = []
    G = []
    B = []
    for i in range(n):
        Ci = colors[z[i]]      
        if type(Ci) == str:
            RGB = CC.to_rgb(Ci)
        else:
            RGB = Ci
        R.append(RGB[0])
        G.append(RGB[1])
        B.append(RGB[2])

    cmap_dict = {}
    cmap_dict['red']   = [(x0[i],R[i],R[i]) for i in range(len(R))]
    cmap_dict['green'] = [(x0[i],G[i],G[i]) for i in range(len(G))]
    cmap_dict['blue']  = [(x0[i],B[i],B[i]) for i in range(len(B))]
    mymap = LinearSegmentedColormap('mymap',cmap_dict)
    return mymap

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make colormap non-linear
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def nonlin_cmap(cmap_in,fac):
    ncol = 256

    orig = pl.get_cmap(cmap_in)(np.linspace(0,1,ncol))
    x0   = np.linspace(0,1,ncol)**fac

    R    = orig[:,0]
    G    = orig[:,1]
    B    = orig[:,2]

    cmap_dict = {}
    cmap_dict['red'] = [(x0[i],R[i],R[i]) for i in range(len(R))]
    cmap_dict['green'] = [(x0[i],G[i],G[i]) for i in range(len(G))]
    cmap_dict['blue'] = [(x0[i],B[i],B[i]) for i in range(len(B))]
    mymap = LinearSegmentedColormap('mymap',cmap_dict)
    return mymap

# -------------------------
# Test colormaps
# -------------------------
def show_cmap(cmap):
    n = 10
    x = np.linspace(0,n,n)
    y = np.linspace(0,n,n)
    z = levs.min() + np.random.random((n,n)) * (levs.max() - levs.min())

    pl.figure()
    pl.contourf(x,y,z,levs,extend='both',cmap=cmap)
    pl.colorbar() 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pre-define some color maps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c1='#7b0000';c2='#e3dab9';c3='#cee1e5';c4='#00037b'
rb      = make_colormap({0:c1,0.499:c2,0.5:c3,1.0:c4})
br      = make_colormap({0:c4,0.499:c3,0.5:c2,1.0:c1})
rwb     = make_colormap({0:c1,0.5:'#ffffff',1.0:c4})
bwr     = make_colormap({0:c4,0.5:'#ffffff',1.0:c1})
red     = make_colormap({0.:'#ffffff',1.0:c1})
blue    = make_colormap({0.:'#ffffff',1.0:c4})
green   = make_colormap({0.:'#ffffff',1.0:'#008600'})
blk     = make_colormap({0.:'#ffffff',1.0:'#000000'})
cld     = make_colormap({0.:'#ffffff',0.2:'#4682b4',0.5:'#0bbd17',0.8:'#fdfe00',1.:'#bf5f2d'})
wnd     = cmap_ncl('precip3_16lev')
rain    = make_colormap({0.:'#ffffff',0.1:'#4682b4',0.2:'#0bbd17',0.4:'#fdfe00',1.:'#bf5f2d'})
wup     = make_colormap({0.:'#ffffff',0.25:'#4682b4',0.5:'#0bbd17',0.75:'#fdfe00',1.:'#cc2900'})
cent    = make_colormap({0.:'#ffffff',0.25:'#fdfe00',0.5:'#0bbd17',0.75:'#fdfe00',1.:'#ffffff'})
wupnl   = make_colormap({0.:'#ffffff',0.1:'#4682b4',0.2:'#0bbd17',0.6:'#fdfe00',1.:'#cc2900'})
rain2   = cmap_ncl('precip2_17lev') 
rain2nl = nonlin_cmap(rain2,2.) 
rain3   = cmap_ncl('precip3_16lev') 
rain3nl = nonlin_cmap(rain3,2.) 
cloud   = pl.cm.Greys

# Test colormaps
if __name__ == "__main__":
    cmaps = ['rb','br','rwb','bwr','red','blue','green','blk','cld','wnd','rain','wup',\
             'cent','wupnl','rain2','rain2nl','rain3','rain3nl','cloud']

    n  = len(cmaps)
    n1 = np.ceil(np.sqrt(n))
    n2 = np.ceil(n / n1)

    print(n,n1,n2)

    data = np.random.random((10,10))

    pl.figure()
    for i in range(n):
        pl.subplot(int(n1),int(n2),i+1)
        pl.imshow(data, interpolation='nearest', cmap=eval(cmaps[i]))
        pl.colorbar()
         
    pl.savefig('color_maps.png')
