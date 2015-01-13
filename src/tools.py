import numpy as np
import pylab as pl
from math import radians, cos, sin, asin, sqrt

R = 6371.   # radius earth

def key_nearest(array,value):
    return (np.abs(array-value)).argmin()

def value_nearest(array,value):
    return (np.abs(array-value)).argmin()

def modplot(ax,minorticks=True,removeax=True,removeaxis=['right','top'],movespine=True,spacing=2):
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
    if('right' in removeaxis):
      ax.spines['right'].set_visible(False)
      ax.get_yaxis().tick_left()
    if('top' in removeaxis):
      ax.spines['top'].set_visible(False)
      ax.get_xaxis().tick_bottom()

"""
Transform coordinate from deg:min:sec to deg.deg
"""
def ctrans(cin):
  temp = cin.split('=')[-1].split(':')    # in case of X=..
  return float(temp[0])+float(temp[1])/60.+float(temp[2])/3600.

"""
Distance (km) between coordinates
"""
def haversine(lon1, lat1, lon2, lat2):
  lon1,lat1,lon2,lat2 = map(radians,[lon1,lat1,lon2,lat2])
  dlon = lon2 - lon1 
  dlat = lat2 - lat1 
  a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
  c = 2 * asin(sqrt(a)) 
  return R * c

"""
Angle from coordinate 1 to 2
"""
def angle(lon1,lat1,lon2,lat2):
  lon1,lat1,lon2,lat2 = map(radians,[lon1,lat1,lon2,lat2])
  dlon = lon2 - lon1 
  dlat = lat2 - lat1 
  y = np.sin(dlon) * np.cos(lat2)
  x = np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(dlon)
  a = np.arctan2(y,x) * 180./np.pi
  return (a+360)%360.

"""
Coordinate given start, distance and direction
"""
def point(lon1,lat1,d,angle):
  angle = angle*np.pi/180.
  lon1,lat1 = map(radians,[lon1,lat1])
  lat2 = np.arcsin(np.sin(lat1)*np.cos(d/R)+np.cos(lat1)*np.sin(d/R)*np.cos(angle))
  lon2 = lon1 + np.arctan2(np.sin(angle)*np.sin(d/R)*np.cos(lat1),np.cos(d/R)-np.sin(lat1)*np.sin(lat2))
  return (lon2*180./np.pi,lat2*180./np.pi)

"""
Simple version distance (km) between coordinates
"""
def distance(lon1,lat1,lon2,lat2):
  lon1,lat1,lon2,lat2 = map(radians,[lon1,lat1,lon2,lat2])
  x = (lon2-lon1)*np.cos((lat1+lat2)/2.)
  y = lat2-lat1
  return sqrt(x*x + y*y) * 6371. 

"""
Simple version angle from coordinate 1 to 2
"""
def angle2(lon1,lat1,lon2,lat2):
  lon1,lat1,lon2,lat2 = map(radians,[lon1,lat1,lon2,lat2])
  a = np.arctan2(lon2-lon1,lat2-lat1) * 180./np.pi 
  return (a+360)%360

