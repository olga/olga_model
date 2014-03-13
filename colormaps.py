import numpy as np
from pylab import *
import glob
import sys
from matplotlib.colors import LinearSegmentedColormap, ColorConverter
from numpy import sort

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create colormap from NCL's .rgb files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def cmap_ncl(name):
  files = glob.glob('colormaps/'+name+'.rgb')
  if(len(files) < 1):
    sys.exit('cant find file %s.*'%name) 

  cols = []
  f = open('colormaps/'+name+'.rgb','r')
  l = 0
  for line in f:
    if(l==0):
      ncolors = int(line.split('=')[-1])
      #print 'found %i colors in %s'%(ncolors,name) 
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create colormap from list of colors
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def make_colormap(colors):
  z = sort(colors.keys())
  n = len(z)
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
  cmap_dict['red'] = [(x0[i],R[i],R[i]) for i in range(len(R))]
  cmap_dict['green'] = [(x0[i],G[i],G[i]) for i in range(len(G))]
  cmap_dict['blue'] = [(x0[i],B[i],B[i]) for i in range(len(B))]
  mymap = LinearSegmentedColormap('mymap',cmap_dict)
  return mymap

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make colormap non-linear
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def nonlin_cmap(cmap_in,fac):
  ncol = 256

  orig = plt.get_cmap(cmap_in)(np.linspace(0,1,ncol))
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pre-define some color maps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c1='#7b0000';c2='#e3dab9';c3='#cee1e5';c4='#00037b'
rb    =  make_colormap({0:c1,0.499:c2,0.5:c3,1.0:c4})
br    =  make_colormap({0:c4,0.499:c3,0.5:c2,1.0:c1})
rwb   = make_colormap({0:c1,0.5:'#ffffff',1.0:c4})
bwr   = make_colormap({0:c4,0.5:'#ffffff',1.0:c1})
red   = make_colormap({0.:'#ffffff',1.0:c1})
blue  = make_colormap({0.:'#ffffff',1.0:c4})
green = make_colormap({0.:'#ffffff',1.0:'#008600'})
blk   = make_colormap({0.:'#ffffff',1.0:'#000000'})
cld   = make_colormap({0.:'#ffffff',0.2:'#4682b4',0.5:'#0bbd17',0.8:'#fdfe00',1.:'#bf5f2d'})
wnd   = make_colormap({0.:'#ffffff',0.4:'#4682b4',0.6:'#0bbd17',0.8:'#fdfe00',1.:'#bf5f2d'})
rain  = make_colormap({0.:'#ffffff',0.1:'#4682b4',0.2:'#0bbd17',0.4:'#fdfe00',1.:'#bf5f2d'})
wup   = make_colormap({0.:'#ffffff',0.25:'#4682b4',0.5:'#0bbd17',0.75:'#fdfe00',1.:'#cc2900'})
cent  = make_colormap({0.:'#ffffff',0.25:'#fdfe00',0.5:'#0bbd17',0.75:'#fdfe00',1.:'#ffffff'})
wupnl = make_colormap({0.:'#ffffff',0.1:'#4682b4',0.2:'#0bbd17',0.6:'#fdfe00',1.:'#cc2900'})
rain2 = cmap_ncl('precip2_17lev') 
rain3 = nonlin_cmap(rain2,2.) 
cloud = nonlin_cmap(cm.Greys_r,2.)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# show colormaps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def show_cmap(cmap_in,nlevs=256):

  n = 30
  c = np.zeros((n,n))
  x = np.linspace(0,6,n)
  y = np.linspace(0,6,n)
  for i in range(n):
    for j in range(n):
      c[i,j] = (np.cos(x[i]) + np.sin(y[j]))/4.+0.5
  
  figure()
  levs = np.linspace(c.min(),c.max(),nlevs)
  contourf(x,y,c,levs,cmap=cmap_in)
  colorbar()



if __name__ == "__main__":
  close('all')
  a=cmap_ncl('WhiteBlueGreenYellowRed')
  show_cmap(a)
  show_cmap(wupnl)
