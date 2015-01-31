import numpy as np
from pylab import *
import urllib2
import os

class readsounding:
    def __init__(self,loc,year,month,day,time):
        # path to sounding (remote)
        path =  'http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT:LIST&YEAR=%04i&MONTH=%02i&FROM=%02i%02i&TO=%02i%02i&STNM=%05i'%(year, month, day, time, day, time, loc)

        try: 
            resp = urllib2.urlopen(path)
            self.success = True
        except: 
            self.success = False
 
        if(self.success): 
            self.p  = []
            self.z  = []
            self.T  = []
            self.Td = []
            self.RH = []
            self.q  = []
            self.V  = []
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
            print('Skipped %2i lines in reading sounding'%skipped)

            self.p  = np.array(self.p )
            self.z  = np.array(self.z )
            self.T  = np.array(self.T )
            self.Td = np.array(self.Td)
            self.RH = np.array(self.RH)
            self.q  = np.array(self.q )
            self.V  = np.array(self.V )
            self.Vd = np.array(self.Vd)

#d = readsounding(10410,2008,12,12,00)
