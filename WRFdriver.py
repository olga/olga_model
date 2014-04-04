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

import numpy as np
import os
import sys
import glob
import urllib
import time
import datetime
import subprocess
import threading
from multiprocessing import Process

#from src.plotdriver import makeplots
from src.settings import *

debug = True

def getGFS(year,month,day,cycle,t0,t1):
    print('Obtaining GFS data...')
    gfsrundir = gfsDataRoot + str(year).zfill(4) + str(month).zfill(2) + str(day).zfill(2) + '/'
    if not os.path.exists(gfsrundir):
        if(debug): print('Making directory %s'%gfsrundir)
        os.mkdir(gfsrundir)
 
    dt = 3   # forecast time interval
 
    nt = int((t1-t0)/dt+1)
    for t in range(nt):
        tact = t0 + t * dt
        loc = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.'+\
                str(year)+\
                str(month).zfill(2)+\
                str(day).zfill(2)+\
                str(cycle).zfill(2)
        fil = 'gfs.t'+\
                str(cycle).zfill(2)+\
                'z.pgrb2f'+\
                str(tact).zfill(2)
        url = loc + '/' + fil

        success = False
        while(success == False):
            if(debug): print('processing %s .. '%fil),
            # Check if file locally available and valid
            if(os.path.isfile(gfsrundir+fil)):
                if(debug): print('found local .. '),
                # Check if same size as remote:
                remote = urllib.urlopen(url)
                meta = remote.info()
                size_remote = meta.getheaders("Content-Length")[0] 
                local = open(gfsrundir+fil,'rb')
                size_local = len(local.read())
                if(int(size_remote) == int(size_local)):
                    if(debug): print('size remote/local match, success!')
                    success = True
                else:
                    if(debug): print('size remote/local differ, re-download .. '),
            # If not, check if available at server:
            if(success == False):
                check = urllib.urlopen(url)
                if(check.code == 200):
                    if(debug): print('file available at GFS server -> downloading')
                    urllib.urlretrieve(url,gfsrundir+fil)
                    # Dont check success, this way the file size check is done again just to be sure
                else:
                    print 'file not found on server, sleep 5min'
                    # File not (yet) available, sleep a while and re-do the checks 
                    time.sleep(300)

    print 'finished GFS at',datetime.datetime.now().time()

def replace(filein,searchstring,value):
    arg =  'sed -i -e "s/\(' +searchstring+ r'\).*/\1 = ' +value+ '/g" ' + filein
    subprocess.call(arg,shell=True)

def updateNamelists(start,end,restart_intv=1440,restart=False,):
    print('updating namelists WPS and WRF...')
    replace(wrfRoot+'namelist.input','start_year' ,str(start.year) +','+str(start.year))
    replace(wrfRoot+'namelist.input','start_month',str(start.month)+','+str(start.month))
    replace(wrfRoot+'namelist.input','start_day'  ,str(start.day)  +','+str(start.day))
    replace(wrfRoot+'namelist.input','start_hour' ,str(start.hour) +','+str(start.hour))
    replace(wrfRoot+'namelist.input','end_year'   ,str(end.year)   +','+str(end.year))
    replace(wrfRoot+'namelist.input','end_month'  ,str(end.month)  +','+str(end.month))
    replace(wrfRoot+'namelist.input','end_day'    ,str(end.day)    +','+str(end.day))
    replace(wrfRoot+'namelist.input','end_hour'   ,str(end.hour)   +','+str(end.hour))

    # Set restart file frequency, and restart flag
    rflag = '.true.' if restart else '.false.'
    replace(wrfRoot+'namelist.input','restart_interval' ,str(restart_intv))
    replace(wrfRoot+'namelist.input',' restart ' ,rflag) # KEEP SPACES'

    replace(wpsRoot+'namelist.wps',  'start_year' ,str(start.year) +','+str(start.year))
    replace(wpsRoot+'namelist.wps',  'start_month',str(start.month)+','+str(start.month))
    replace(wpsRoot+'namelist.wps',  'start_day'  ,str(start.day)  +','+str(start.day))
    replace(wpsRoot+'namelist.wps',  'start_hour' ,str(start.hour) +','+str(start.hour))
    replace(wpsRoot+'namelist.wps',  'end_year'   ,str(end.year)   +','+str(end.year))
    replace(wpsRoot+'namelist.wps',  'end_month'  ,str(end.month)  +','+str(end.month))
    replace(wpsRoot+'namelist.wps',  'end_day'    ,str(end.day)    +','+str(end.day))
    replace(wpsRoot+'namelist.wps',  'end_hour'   ,str(end.hour)   +','+str(end.hour))

def runWPS(start,clean=True):
    print('Running WPS...')
    # Grr, we have to call the routines from the directory itself...
    os.chdir(wpsRoot)

    if(clean):
        # Cleanup stuff from previous day
        subprocess.call('rm GRIBFILE*',shell=True)
        subprocess.call('rm FILE*',shell=True)
        subprocess.call('rm met_em*',shell=True)

    # Each log will be saved in olgaLogs with yyyymmddhh added
    logappend = str(start.year).zfill(4)  +\
                str(start.month).zfill(2) +\
                str(start.day).zfill(2)   +\
                str(start.hour).zfill(2)

    # Make directory for the logs, if it doesn't exist
    if not os.path.exists(olgaLogs):
        os.mkdir(olgaLogs)

    # Run geogrid
    if(debug): print('... WPS -> geogrid')
    subprocess.call('./geogrid.exe >& %sgeogrid.%s'%(olgaLogs,logappend),shell=True)

    # Link the GFS data to the WPS directory
    gfsData = gfsDataRoot + str(start.year).zfill(4) + str(start.month).zfill(2) + str(start.day).zfill(2)
    subprocess.call('./link_grib.csh '+gfsData+'/gfs*',shell=True)

    # Check if Vtable present, if not link it -> unpack GFS data
    if(debug): print('... WPS -> ungrib')
    if(not os.path.isfile('Vtable')):
        subprocess.call('ln -s ungrib/Variable_Tables/Vtable.GFS Vtable',shell=True)
    subprocess.call('./ungrib.exe >& %sungrib.%s'%(olgaLogs,logappend),shell=True)

    # Interpolate to correct grid
    if(debug): print('... WPS -> metgrid')
    subprocess.call('./metgrid.exe >& %smetgrid.%s'%(olgaLogs,logappend),shell=True)

    os.chdir(olgaRoot)

def runWRF(start,clean=False):
    print('Running WRF...')
    # Grr, we have to call the routines from the directory itself...
    os.chdir(wrfRoot)

    # Remove restart file from previous days
    if(clean):
        subprocess.call('rm wrfrst*',shell=True)

    # Each log will be saved in olgaLogs with yyyymmddhh added
    logappend = str(start.year).zfill(4)  +\
                str(start.month).zfill(2) +\
                str(start.day).zfill(2)   +\
                str(start.hour).zfill(2)

    # Make directory for the logs, if it doesn't exist
    if not os.path.exists(olgaLogs):
        os.mkdir(olgaLogs)
 
    # Link the met_em input files 
    subprocess.call('rm met_em*',shell=True)
    subprocess.call('ln -s '+wpsRoot+'met_em* .',shell=True)
    # Run real
    if(debug): print('... WRF -> real.exe')
    subprocess.call('./real.exe >& %sreal.%s'%(olgaLogs,logappend),shell=True)
    # Run WRF as background process to allow postprocessing to run at the same time..
    if(debug): print('... WRF -> wrf.exe')
    subprocess.call('./wrf.exe >& %swrf.%s &'%(olgaLogs,logappend),shell=True)

    os.chdir(olgaRoot)

"""
Wait untill the correct restart file is available
"""
def waitforWRF(start,end):
    print('Waiting for WRF to finish')
    wrfrst = wrfRoot + 'wrfrst_d01_' + str(end.year).zfill(4)  + \
                                 '-' + str(end.month).zfill(2) + \
                                 '-' + str(end.day).zfill(2) + '_00:00:00'
  
    while(True):
        if(not os.path.isfile(wrfrst)):
            time.sleep(5)
        else: 
            print 'finished WRF at',datetime.datetime.now().time()
            break

"""
asdfasdf
"""
def moveWRFout(start,basestr,nrun):
    wrfouts = glob.glob(wrfRoot+'wrfout_d0*')

    if(len(wrfouts)>0):
        for wrfout in wrfouts:
            domain = int(wrfout.split('_')[1][-1])
            outname = '%s_d%i_%i.nc'%(basestr,domain,nrun)
            subprocess.call('mv %s %s%s'%(wrfout,wrfDataRoot,outname),shell=True)

"""
asdfasdf
"""
def triggerPlots(start,cycle,nrun):
    date = str(start.year).zfill(4) + str(start.month).zfill(2) + str(start.day).zfill(2)

    subprocess.call('python2 %s/src/plotdriver.py %s %i %i %i %i 2 time >& /dev/null &'%\
                   (olgaRoot,start.year,start.month,start.day,nrun,cycle),shell=True)
    subprocess.call('python2 %s/src/plotdriver.py %s %i %i %i %i 2 sounding >& /dev/null &'%\
                   (olgaRoot,start.year,start.month,start.day,nrun,cycle),shell=True)
    subprocess.call('python2 %s/src/plotdriver.py %s %i %i %i %i 1 maps >& /dev/null  &'%\
                   (olgaRoot,start.year,start.month,start.day,nrun,cycle),shell=True)
    subprocess.call('python2 %s/src/plotdriver.py %s %i %i %i %i 2 maps >& /dev/null '%\
                   (olgaRoot,start.year,start.month,start.day,nrun,cycle),shell=True)

    print 'finished plots at',datetime.datetime.now().time()


if __name__ == "__main__":
    
    # Get command line arguments
    modes = ['all','post']
    if(len(sys.argv) != 2):
        sys.exit('provide mode: {all,post}')
    else:
        mode = sys.argv[1]
        if(mode not in modes):
            sys.exit('mode %s invalid'%mode)

    # Get current date   
    year     = '2014' #time.strftime('%Y')
    month    = '04' #time.strftime('%m')
    day      = '01' #time.strftime('%d')
    tstart   = 00   #  

    cycle    = 0    # which GFS cycle? {0,6,12,18}
    dtinput  = 3    # input dt of GFS (==3)
    ttotal   = 48   # total number of hours to simulate
    tfirst   = 24   # initial hours to simulate
    ndom     = 2    # number of domains

    basestr  = "%s%s%s_t%02iz"%(year,month,day,cycle)

    startstruct  = datetime.datetime.strptime('%s %s %s %s'%(day,month,year,tstart),"%d %m %Y %H")
    endstruct1   = startstruct + datetime.timedelta(hours=tfirst)
    endstruct2   = startstruct + datetime.timedelta(hours=ttotal)

    if(mode=='all'):
        # Start first 'tfirst' hours of sequence
        getGFS(year,month,day,cycle=cycle,t0=tstart,t1=tfirst)
        updateNamelists(startstruct,endstruct1,tfirst*60,False)
        runWPS(startstruct,clean=True)
        runWRF(startstruct,clean=True)

        # While WRF is running, download rest of GFS
        getGFS(year,month,day,cycle=cycle,t0=tfirst,t1=ttotal)

        # Wait untill the restart file is available
        waitforWRF(startstruct,endstruct1)

    if(mode=='all' or mode =='post'):
        moveWRFout(startstruct,basestr,nrun=1)
        triggerPlots(startstruct,cycle,nrun=1)

    if(mode=='all'):
        # Run remaining part of simulation
        updateNamelists(endstruct1,endstruct2,9999,True)
        runWPS(startstruct,clean=False)
        runWRF(startstruct,clean=False)
        waitforWRF(startstruct,endstruct2)

    if(mode=='all' or mode =='post'):
        moveWRFout(startstruct,basestr,nrun=2)
        triggerPlots(startstruct,cycle,nrun=2)



