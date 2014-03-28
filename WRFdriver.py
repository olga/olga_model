import numpy as np
import os
import urllib
import time
import subprocess
import threading
from netCDF4 import Dataset

# -------------------
# User settings
# -------------------
olgaRoot = '/home/bart/meteo/WRFnl/olga/'
gfsDataRoot = '/home/bart/meteo/WRFnl/GFSdata/'
wpsRoot = '/home/bart/meteo/WRFnl/WPSv531/'
wrfRoot = '/home/bart/meteo/WRFnl/WRFv351/run/'

debug = True

def getGFS(year,month,day,cycle,t0,t1):
    print('Obtaining GFS data...')
    dirout = gfsDataRoot + str(year).zfill(4) + str(month).zfill(2) + str(day).zfill(2) + '/'
    if not os.path.exists(gfsDataRoot):
        if(debug): print('Making directory %s'%gfsDataRoot)
        os.mkdir(dirout)
 
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

        if(debug): print('processing file %s'%fil)

        success = False
        while(success == False):
            # Check if file locally available and valid
            if(os.path.isfile(gfsDataRoot+fil)):
                if(debug): print('found %s on local system'%fil)
                # Check if same size as remote:
                remote = urllib.urlopen(url)
                meta = remote.info()
                size_remote = meta.getheaders("Content-Length")[0] 
                local = open(gfsDataRoot+fil,'rb')
                size_local = len(local.read())
                if(int(size_remote) == int(size_local)):
                    if(debug): print('file size remote and local match, success!')
                    success = True
                else:
                    if(debug): print('file size remote and local differ, re-download')
            # If not, check if available at server:
            if(success == False):
                check = urllib.urlopen(url)
                if(check.code == 200):
                    if(debug): print('file available at GFS server -> downloading')
                    urllib.urlretrieve(url,gfsDataRoot+fil)
                else:
                    # Wait a while and re-do checks
                    time.sleep(300)

def replace(filein,searchstring,value):
    arg =  'sed -i -e "s/\(' +searchstring+ r'\).*/\1 = ' +value+ '/g" ' + filein
    subprocess.call(arg,shell=True)

def updateNamelists(year,month,day,dt):
    print('updating namelists WPS and WRF...')
    # 1. update namelist.wps and namelist.input
    replace(wrfRoot+'namelist.input','start_year',str(year)+','+str(year))
    replace(wrfRoot+'namelist.input','start_month',str(month)+','+str(month))
    replace(wrfRoot+'namelist.input','start_day',str(day)+','+str(day))
    replace(wrfRoot+'namelist.input','end_year',str(year)+','+str(year))
    replace(wrfRoot+'namelist.input','end_month',str(month)+','+str(month))
    replace(wrfRoot+'namelist.input','end_day',str(int(day)+dt)+','+str(int(day)+dt))

    replace(wpsRoot+'namelist.wps','start_year',str(year)+','+str(year))
    replace(wpsRoot+'namelist.wps','start_month',str(month)+','+str(month))
    replace(wpsRoot+'namelist.wps','start_day',str(day)+','+str(day))
    replace(wpsRoot+'namelist.wps','end_year',str(year)+','+str(year))
    replace(wpsRoot+'namelist.wps','end_month',str(month)+','+str(month))
    replace(wpsRoot+'namelist.wps','end_day',str(int(day)+dt)+','+str(int(day)+dt))

def runWPS():
    print('Running WPS...')
    # Grr, we have to call the routines from the directory itself...
    os.chdir(wpsRoot)

    # Cleanup stuff from previous day
    subprocess.call('rm GRIBFILE*',shell=True)
    subprocess.call('rm FILE*',shell=True)
    subprocess.call('rm met_em*',shell=True)

    # Run geogrid
    if(debug): print('... WPS -> geogrid')
    subprocess.call('./geogrid.exe >& geogrid.log',shell=True)
    subprocess.call('./link_grib.csh '+gfsDataRoot+'gfs*',shell=True)
    if(debug): print('... WPS -> ungrib')
    # Check if Vtable present, if not link it
    if(not os.path.isfile('Vtable')):
        subprocess.call('ln -s ungrib/Variable_Tables/Vtable.GFS Vtable',shell=True)
    subprocess.call('./ungrib.exe >& ungrib.log',shell=True)
    if(debug): print('... WPS -> metgrid')
    subprocess.call('./metgrid.exe >& metgrid.log',shell=True)

    os.chdir(olgaRoot)

def runWRF():
    print('Running WRF...')
    # Grr, we have to call the routines from the directory itself...
    os.chdir(wrfRoot)
 
    # Link the met_em input files 
    subprocess.call('rm met_em*',shell=True)
    subprocess.call('ln -s '+wpsRoot+'met_em* .',shell=True)
    # Run real
    if(debug): print('... WRF -> real.exe')
    subprocess.call('./real.exe >& real.log',shell=True)
    # Run WRF as backgroudn process to allow postprocessing to run at the same time..
    if(debug): print('... WRF -> wrf.exe')
    subprocess.call('./wrf.exe >& wrf.log &',shell=True)

    os.chdir(olgaRoot)

def wait4WRF(year,month,day):
    print('Waiting for WRF to finish')
    wrfout = wrfRoot + 'wrfout_d01_' + str(year).zfill(4) + \
                                 '-' + str(month).zfill(2) + \
                                 '-' + str(day).zfill(2) + '_00:00:00'

    if(os.path.isfile(wrfout)):
        nc = Dataset(wrfout,'r')
        time = nc.variables["Times"][:]
        print(time)

if __name__ == "__main__":
    # Get current date   
    year = 2014 #time.strftime('%Y')
    month = 03 #time.strftime('%m')
    day = 15 #time.strftime('%d')

    #getGFS(year,month,day,cycle=0,t0=0,t1=24)
    #updateNamelists(year,month,day,1)
    #runWPS()
    #runWRF()
    wait4WRF(year,month,day)



    #postproc()


    # Download GFS data
    #t1a = threading.Thread(target=getGFS(year,month,day,cycle=0,t0=0,t1=24))
    # Update namelists for first day of forecast
    #t2a = threading.Thread(target=updateNamelists(year,month,day,1))

    #t1a.start()  # Download first day GFS
    #t2a.start()  # Update namelist for 1st stage
    #t1a.join()   # Wait for download to finish
    #t2a.join()   # Wait for updateing namelist to finish

    #runWPS()
    #runWRF()

    ## Start WPS, and start downloading rest of GFS
    #t3 = threading.Thread(target=runWPS)
    #t3.start()
    #t1b.start()
    #t1c.start()
   
    #  
    #t3.join()
    
     
    






