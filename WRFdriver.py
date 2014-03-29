import numpy as np
import os
import urllib
import time
import datetime
import subprocess
import threading
from multiprocessing import Process
from netCDF4 import Dataset

# -------------------
# User settings
# -------------------
olgaRoot    = '/home/bart/meteo/WRFnl/olga/' # Full path to OLGA scripts
olgaLogs    = olgaRoot+'/logs/' # Location to save logs
gfsDataRoot = '/home/bart/meteo/WRFnl/GFSdata/' # Path to root of GFS data
wpsRoot     = '/home/bart/meteo/WRFnl/WPSv531/' # Path to root of WPS
wrfRoot     = '/home/bart/meteo/WRFnl/WRFv351/run/' # Path to root of WRF

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
                    # File not (yet) available, sleep a while and re-do the checks 
                    time.sleep(300)

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
 
    # Link the met_em input files 
    subprocess.call('rm met_em*',shell=True)
    subprocess.call('ln -s '+wpsRoot+'met_em* .',shell=True)
    # Run real
    if(debug): print('... WRF -> real.exe')
    subprocess.call('./real.exe >& %sreal.%s'%(olgaLogs,logappend),shell=True)
    # Run WRF as backgroudn process to allow postprocessing to run at the same time..
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
            break

if __name__ == "__main__":
    # Get current date   
    year     = 2014 #time.strftime('%Y')
    month    = 03   #time.strftime('%m')
    day      = 15   #time.strftime('%d')
    tstart   = 00   #  

    cycle    = 0    # which GFS cycle? {0,6,12,18}
    dtinput  = 3    # input dt of GFS (==3)
    ttotal   = 48   # total number of hours to simulate
    tfirst   = 24   # initial hours to simulate

    startstruct  = datetime.datetime.strptime('%i %i %i %i'%(day,month,year,tstart),"%d %m %Y %H")
    endstruct1   = startstruct + datetime.timedelta(hours=tfirst)
    endstruct2   = startstruct + datetime.timedelta(hours=ttotal)

    # Start first 'tfirst' hours of sequence
    getGFS(year,month,day,cycle=cycle,t0=tstart,t1=tfirst)
    updateNamelists(startstruct,endstruct1,tfirst*60,False)
    runWPS(startstruct,clean=True)
    runWRF(startstruct,clean=True)

    # While WRF is running, download rest of GFS
    getGFS(year,month,day,cycle=cycle,t0=tfirst,t1=ttotal)

    # Wait untill the restart file is available
    waitforWRF(startstruct,endstruct1)

    # Run remaining part of simulation
    updateNamelists(endstruct1,endstruct2,9999,True)
    runWPS(startstruct,clean=False)
    runWRF(startstruct,clean=False)


    """
    Goal was to run certain processes (downloading, WRF) parallel.. Doesnt work?
    """
    #t1a = Process(target=getGFS(year,month,day,cycle=cycle,t0=tstart,t1=tfirst))
    #t1a.start() # Start downloading first day of GFS 
    #t2a = Process(target=updateNamelists(startstruct,endstruct1,tfirst*60,False))
    #t2a.start() # Update namelists

    #t1a.join() # Wait for download to finish
    #t2a.join() # Wait for updating namelists to finish

    #t3a = Process(target=runWPS(startstruct,clean=True))
    #t3a.start() # Start WPS for first day
    #t1b = Process(target=getGFS(year,month,day,cycle=cycle,t0=tfirst,t1=ttotal))
    #t1b.start() # Continue downloading GFS

    #t3a.join() # Wait for WPS to finish
    #t4a = Process(target=runWRF(startstruct,clean=True))
    #t4a.start() # Start WRF for first day
    #t5a = Process(target=waitforWRF(startstruct,endstruct1))
    #t5a.start() # Start waiting for WRF to finish

    #t4a.join() # WRF return directly....
    #t5a.join() # End of first part of simulation

    #### POSTPROC

    #t2b = Process(target=updateNamelists(endstruct1,endstruct2,9999,True))
    #t2b.start() # Update namelists
    #t1b.join() # GFS needs to be finished
    #t2b.join() # Namelists finished
    #
    #t3b = Process(target=runWPS(startstruct,clean=False))
    #t3b.start()
    #t3b.join()
 
    #t4b = Process(target=runWRF(startstruct,clean=False))
    #t4b.start()
    #t5b = Process(target=waitforWRF(startstruct,endstruct2))
    #t5b.start()
    #t4b.join()
    #t5b.join()

    ### POSTPROC

    
     
    






