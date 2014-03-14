import numpy as np
from pylab import *
import urllib2
import os

class readsounding:
  def __init__(self,name,loc,year,month,day,time):

    self.name = name
    loc = str(loc).zfill(2)
    year = str(year)
    month = str(month).zfill(2)
    day = str(day).zfill(2)
    time = str(time).zfill(2)
    
    # path to sounding (remote)
    path =  'http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&YEAR='+year+'&MONTH='+month+'&FROM='+day+time+'&TO='+day+time+'&STNM='+loc

    try: 
      resp = urllib2.urlopen(path)
    except: 
      sys.exit('file not local, cant retrieve from uwyo.edu')
 
    self.p = []
    self.z = []
    self.T = []
    self.Td = []
    self.RH = []
    self.q = []
    self.V = []
    self.Vd = []

    # Holy @#$@$# how ugly
    sound = resp.read().split('\n')
    skipped = 0
    for i in range(10,len(sound)):
      if(len(sound[i].split())==11):
        try:
          vals = sound[i].split()
          self.p.append(float(vals[0]))
          self.z.append(float(vals[1]))
          self.T.append(float(vals[2]))
          self.Td.append(float(vals[3]))
          self.RH.append(float(vals[4]))
          self.q.append(float(vals[5]))
          self.V.append(float(vals[7]))
          self.Vd.append(float(vals[6]))
        except ValueError:
          skipped += 1
    print 'Skipped %2i lines in reading sounding'%skipped

    self.p  = array( self.p  )
    self.z  = array( self.z  )
    self.T  = array( self.T  )
    self.Td = array( self.Td )
    self.RH = array( self.RH )
    self.q  = array( self.q  )
    self.V  = array( self.V  )
    self.Vd = array( self.Vd )

#d = readsounding(10410,2008,12,12,00)
