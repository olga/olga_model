import numpy as np
from pylab import *

def key_nearest(array,value):
    return (np.abs(array-value)).argmin()

def value_nearest(array,value):
    return (np.abs(array-value)).argmin()

def modplot(ax,minorticks=True,removeax=True,movespine=True,spacing=2):
  if(minorticks):
    from matplotlib.ticker import AutoMinorLocator
    minorLocator   = AutoMinorLocator()
    ax.yaxis.set_minor_locator(minorLocator)
    minorLocator   = AutoMinorLocator()
    ax.xaxis.set_minor_locator(minorLocator)

  if(movespine):
    for loc, spine in ax.spines.items():
      spine.set_position(('outward',spacing)) # outward by 10 points

  if(removeax):
    ax.spines['right'].set_visible(False)
    ax.get_yaxis().tick_left()
    ax.spines['top'].set_visible(False)
    ax.get_xaxis().tick_bottom()

def cmap_discrete(cmap_in,ints):
  return plt.get_cmap(cmap_in)(ints)


