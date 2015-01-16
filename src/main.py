
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

from src.plotdriver import plot_driver

debug = True  

## Downloads the requested GFS data, or waits until available
# @param olga Pointer to object with OLGA settings
def downloadGFS(olga,islice):
    print('Obtaining GFS data...')
    gfsbase = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs'
    gfsrundir = '%s%04i%02i%02i/'%(olga.gfsDataRoot,olga.year,olga.month,olga.day) # where to save files

    dtGFS       = 3. # Time step of GFS input data [h]

    # Check if GFS data directory exists, if not create it
    if not os.path.exists(gfsrundir):
        if(debug): print('Making directory %s'%gfsrundir)
        os.mkdir(gfsrundir)

    # Calculate first and last hour to download
    t0 = olga.tstart + (islice+0) * olga.tslice
    t1 = olga.tstart + (islice+1) * olga.tslice
    nt = int((t1-t0) / dtGFS) + 1 

    # Loop over time steps
    for t in range(nt):
        tact = int(t0 + t * dtGFS) # Forecast time
        loc = '%s.%04i%02i%02i%02i'%(gfsbase,olga.year,olga.month,olga.day,olga.cycle) # Location at server
        fil = 'gfs.t%02iz.pgrb2f%02i'%(olga.cycle,tact) # File at server
        url = '%s/%s'%(loc,fil) # Path to file at server

        success = False
        while(success == False):
            # Test: put try around everything to catch weird exceptions
            try:
                if(debug): print('processing %s .. '%fil),
                # Check if file locally available and valid
                if(os.path.isfile(gfsrundir+fil)):
                    if(debug): print('found local .. '),
                    # Check if same size as remote:
                    remote = urllib.urlopen(url)
                    # for re-running old cases, the file might have been removed online.
                    if(remote.code == 200):
                        meta = remote.info()
                        size_remote = meta.getheaders("Content-Length")[0] 
                        local = open(gfsrundir+fil,'rb')
                        size_local = len(local.read())
                        if(int(size_remote) == int(size_local)):
                            if(debug): print('size remote/local match, success!')
                            success = True
                        else:
                            if(debug): print('size remote/local differ, re-download .. '),
                    else:
                        if(debug): print('remote not available for comparision: assume all is fine :)')
                        success = True
                # If not, check if available at server:
                if(success == False):
                    check = urllib.urlopen(url)
                    if(check.code == 200):
                        # File available, download! 
                        if(debug): print('file available at GFS server -> downloading')
                        urllib.urlretrieve(url,gfsrundir+fil)
                    else:
                        # File not (yet) available, sleep a while and re-do the checks 
                        print('file not found on server, sleep 5min')
                        time.sleep(300)
            except:
                # Something weird happened. Sleep 5 minutes, try again
                print('weird exception:')
                print(sys.exc_info()[0]) 
                time.sleep(300)

    print('finished GFS at %s'%datetime.datetime.now().time())

## replace everything after '=' on line containing 'searchstring' with 'value' in 'filein' 
# @param filein Path to file
# @param searchstring String to search for
# @param value Value to replace
def replace(filein,searchstring,value):
    arg =  'sed -i -e "s/\(' +searchstring+ r'\).*/\1 = ' +value+ '/g" ' + filein
    subprocess.call(arg,shell=True,executable='/bin/bash')

## Slow way of creating a string consisting of n times the same string...
# @param string String to glue
# @param n How many times to glue
# @param separator Separator in string
def printn(string,n,separator=','):
    str_out = ''
    for i in range(n):
        str_out += str(string) + str(separator)
    return str_out

## Update WRF and WPS namelists
# @param olga Pointer to object with OLGA settings
def updateNamelists(olga):
    print('updating namelists WPS and WRF...')

    # Update WRF namelist
    replace(olga.wrfRoot+'namelist.input','start_year',   printn(olga.startstruct.year,   olga.ndom))
    replace(olga.wrfRoot+'namelist.input','start_month',  printn(olga.startstruct.month,  olga.ndom))
    replace(olga.wrfRoot+'namelist.input','start_day',    printn(olga.startstruct.day,    olga.ndom))
    replace(olga.wrfRoot+'namelist.input','start_hour',   printn(olga.startstruct.hour,   olga.ndom))
    replace(olga.wrfRoot+'namelist.input','start_minute', printn(olga.startstruct.minute, olga.ndom))
    replace(olga.wrfRoot+'namelist.input','start_second', printn(olga.startstruct.second, olga.ndom))
    replace(olga.wrfRoot+'namelist.input','end_year',     printn(olga.endstruct.year,     olga.ndom))
    replace(olga.wrfRoot+'namelist.input','end_month',    printn(olga.endstruct.month,    olga.ndom))
    replace(olga.wrfRoot+'namelist.input','end_day',      printn(olga.endstruct.day,      olga.ndom))
    replace(olga.wrfRoot+'namelist.input','end_hour',     printn(olga.endstruct.hour,     olga.ndom))
    replace(olga.wrfRoot+'namelist.input','end_minute',   printn(olga.endstruct.minute,   olga.ndom))
    replace(olga.wrfRoot+'namelist.input','end_second',   printn(olga.endstruct.second,   olga.ndom))

    # Set restart file frequency, restart flag and number of domains
    rflag = '.true.' if olga.islice>0 else '.false.'
    replace(olga.wrfRoot+'namelist.input','restart_interval' ,str(olga.tslice*60))
    replace(olga.wrfRoot+'namelist.input',' restart ',   rflag) # KEEP SPACES AROUND restart'
    replace(olga.wrfRoot+'namelist.input','max_dom',     str(olga.ndom))

    # Update WPS namelist
    replace(olga.wpsRoot+'namelist.wps','start_year',    printn(olga.startstruct.year,   olga.ndom))
    replace(olga.wpsRoot+'namelist.wps','start_month',   printn(olga.startstruct.month,  olga.ndom))
    replace(olga.wpsRoot+'namelist.wps','start_day',     printn(olga.startstruct.day,    olga.ndom))
    replace(olga.wpsRoot+'namelist.wps','start_hour',    printn(olga.startstruct.hour,   olga.ndom))
    replace(olga.wpsRoot+'namelist.wps','start_minute',  printn(olga.startstruct.minute, olga.ndom))
    replace(olga.wpsRoot+'namelist.wps','start_second',  printn(olga.startstruct.second, olga.ndom))
    replace(olga.wpsRoot+'namelist.wps','end_year',      printn(olga.endstruct.year,     olga.ndom))
    replace(olga.wpsRoot+'namelist.wps','end_month',     printn(olga.endstruct.month,    olga.ndom))
    replace(olga.wpsRoot+'namelist.wps','end_day',       printn(olga.endstruct.day,      olga.ndom))
    replace(olga.wpsRoot+'namelist.wps','end_hour',      printn(olga.endstruct.hour,     olga.ndom))
    replace(olga.wpsRoot+'namelist.wps','end_minute',    printn(olga.endstruct.minute,   olga.ndom))
    replace(olga.wpsRoot+'namelist.wps','end_second',    printn(olga.endstruct.second,   olga.ndom))
    replace(olga.wpsRoot+'namelist.wps','max_dom',       str(olga.ndom))

## Run the WPS steps
# @param olga Pointer to object with OLGA settings
def execWPS(olga):
    print('Running WPS at %s'%datetime.datetime.now().time())
    # Grr, we have to call the routines from the directory itself...
    os.chdir(olga.wpsRoot)

    if(olga.islice==0):
        # Cleanup stuff from previous day
        subprocess.call('rm GRIBFILE* >& /dev/null',shell=True,executable='/bin/bash')
        subprocess.call('rm FILE* >& /dev/null',shell=True,executable='/bin/bash')
        subprocess.call('rm met_em* >& /dev/null',shell=True,executable='/bin/bash')

    # Each log will be saved in olgaLogs with yyyymmddhh added
    ss = olga.startstruct
    logappend = '%04i%02i%02i%02i_r%i'%(ss.year,ss.month,ss.day,ss.hour,olga.islice)

    # Make directory for the logs, if it doesn't exist
    if not os.path.exists(olga.olgaLogs):
        os.mkdir(olga.olgaLogs)

    # Run geogrid
    if(debug): print('... WPS -> geogrid')
    subprocess.call('./geogrid.exe >& %sgeogrid.%s'%(olga.olgaLogs,logappend),shell=True,executable='/bin/bash')

    # Link the GFS data to the WPS directory
    gfsData = '%s%04i%02i%02i'%(olga.gfsDataRoot,olga.year,olga.month,olga.day)
    subprocess.call('./link_grib.csh '+gfsData+'/gfs*',shell=True,executable='/bin/bash')

    # Check if Vtable present, if not link it
    if(debug): print('... WPS -> ungrib')

    # Run ungrib
    subprocess.call('./ungrib.exe >& %sungrib.%s'%(olga.olgaLogs,logappend),shell=True,executable='/bin/bash')

    # Run metgrid
    if(debug): print('... WPS -> metgrid')
    subprocess.call('./metgrid.exe >& %smetgrid.%s'%(olga.olgaLogs,logappend),shell=True,executable='/bin/bash')

    print('finished WPS at %s'%datetime.datetime.now().time())
    os.chdir(olga.domainRoot)

## Run the WRF steps
# @param olga Pointer to object with OLGA settings
def execWRF(olga):
    print('Running WRF at %s'%datetime.datetime.now().time())
    # Grr, we have to call the routines from the directory itself...
    os.chdir(olga.wrfRoot)

    # Remove restart file and output from previous day(s)
    if(olga.islice==0):
        subprocess.call('rm wrfrst* >& /dev/null',shell=True,executable='/bin/bash')
        subprocess.call('rm wrfout* >& /dev/null',shell=True,executable='/bin/bash')

    # Each log will be saved in olgaLogs with yyyymmddhh added
    ss = olga.startstruct
    logappend = '%04i%02i%02i%02i_r%i'%(ss.year,ss.month,ss.day,ss.hour,olga.islice)

    # Make directory for the logs, if it doesn't exist
    if not os.path.exists(olga.olgaLogs):
        os.mkdir(olga.olgaLogs)
 
    # Link the met_em input files 
    subprocess.call('rm met_em*',shell=True,executable='/bin/bash')
    subprocess.call('ln -s '+olga.wpsRoot+'met_em* .',shell=True,executable='/bin/bash')

    if(olga.ompThreads > 1):
        subprocess.call('export OMP_NUM_THREADS=%i'%(olga.omp_thr),shell=True,executable='/bin/bash')

    # Run real
    #if(debug): print('... WRF -> real.exe')
    #if(olga.mpiTasks > 1):
    #    subprocess.call('mpirun -n %i ./real.exe >& %sreal.%s'%(olga.mpiTasks,olga.olgaLogs,logappend),shell=True,executable='/bin/bash')
    #else:
    #    subprocess.call('./real.exe >& %sreal.%s'%(olga.olgaLogs,logappend),shell=True,executable='/bin/bash')

    # Run WRF as background process to allow other processes (download GFS, ..) to run at the same time..
    #if(debug): print('... WRF -> wrf.exe')
    #if(olga.mpiTasks > 1):
    #    subprocess.call('mpirun -n %i ./wrf.exe >& %swrf.%s &'%(olga.mpiTasks,olga.olgaLogs,logappend),shell=True,executable='/bin/bash')
    #else:
    #    subprocess.call('./wrf.exe >& %swrf.%s &'%(olga.olgaLogs,logappend),shell=True,executable='/bin/bash')

    # Load balancer Thunder
    subprocess.call('sbatch run.slurm',shell=True,executable='/bin/bash')

    os.chdir(olga.domainRoot)

## Wait until the required restart file is available (i.e. WRF finished)
# @param olga Pointer to object with OLGA settings
def wait4WRF(olga):
    print('Waiting for WRF to finish')

    es = olga.endstruct
    wrfrst = '%swrfrst_d01_%04i-%02i-%02i_%02i:%02i:00'%(olga.wrfRoot,es.year,es.month,es.day,es.hour,es.minute)

    # Check if 'wrfrst' is available, else sleep
    while(True):
        if(not os.path.isfile(wrfrst)):
            time.sleep(10)
        else: 
            print('finished WRF at %s'%datetime.datetime.now().time())
            break

## Copy WRF output to wrfDataRoot
# @param olga Pointer to object with OLGA settings
def moveWRFOutput(olga):
    ss = olga.startstruct
    for dom in range(olga.ndom):
        outname = '%04i%02i%02i_t%02iz_d%i.nc'%(olga.year,olga.month,olga.day,olga.cycle,dom+1)
        wrfouts = glob.glob('%swrfout_d0%i*'%(olga.wrfRoot,dom+1))

        if(np.size(wrfouts) > 0):
            wrfouts.sort() # sort to make sure that the time order is correct in the merged file
            tmp = ''
            for i in range(np.size(wrfouts)):
                tmp += wrfouts[i] + ' '

            # For now, use ncrcat from NCO to merge the files.
            # If this turns out to be problematic (availability NCO on different linux distributions),
            # write own routine to do the merge (shouldn't be difficult)
            subprocess.call('ncrcat -O %s %s%s'%(tmp,olga.wrfDataRoot,outname),shell=True,executable='/bin/bash')

## Create the plots / maps / soundings
# @param olga Pointer to object with OLGA settings
def execPlots(olga):
    print('starting plots at %s'%datetime.datetime.now().time())

    # Spawn different processes for each mape type, saves quite some time..
    for dom in range(olga.ndom):
        if(olga.maps[dom]==True):
            print('making maps for domain=%i'%(dom+1))
            #plot_driver(olga,dom,'maps')
            pmap = Process(target=plot_driver, args=(olga,dom,'maps',))
            pmap.start()

        if(olga.meteogr[dom]):
            print('making time series for domain=%i'%(dom+1))
            #plot_driver(olga,dom,'time')
            pmgram = Process(target=plot_driver, args=(olga,dom,'time',))
            pmgram.start()

        if(olga.sounding[dom]==True):
            print('making soundings for domain=%i'%(dom+1))
            #plot_driver(olga,dom,'sounding')
            psound = Process(target=plot_driver, args=(olga,dom,'sounding',))
            psound.start()

        if(olga.maps[dom]==True):
            pmap.join()
        if(olga.meteogr[dom]):
            pmgram.join()
        if(olga.sounding[dom]==True):
            psound.join()

    print('finished plots at %s'%datetime.datetime.now().time())
