
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
res = '0.25' # '0.25' or 'old'

## Execute command
# @param task Task to start
def execute(task):
    subprocess.call(task, shell=True, executable='/bin/bash')

def progress(count, blockSize, totalSize):
    percent = int(count*blockSize*100/totalSize)
    sys.stdout.write("\r ... %d percent"%percent)
    sys.stdout.flush()

## Print to stdout, and flush buffer
def printf(message):
    print(message)
    sys.stdout.flush() 

## Download single GFS input file, or waits until available
# @ param blabla
def downloadGFSFile(localfile, remoteurl, filename):
    if(debug): printf('processing %s .. '%filename)
    success = False
    while(success == False):
        # Test: put try around everything to catch weird exceptions
        try:
            # Check if file locally available and valid
            if(os.path.isfile(localfile)):
                if(debug): printf('found %s local'%filename)
                # Check if same size as remote:
                remote = urllib.urlopen(remoteurl)
                # for re-running old cases, the file might have been removed online.
                if(remote.code == 200):
                    meta = remote.info()
                    size_remote = meta.getheaders("Content-Length")[0] 
                    local = open(localfile,'rb')
                    size_local = len(local.read())
                    if(int(size_remote) == int(size_local)):
                        if(debug): printf('size %s remote/local match, success!'%filename)
                        success = True
                    else:
                        if(debug): printf('size %s remote/local differ, re-download .. '%filename)
                else:
                    if(debug): printf('remote %s not available (probably old run), pass'%filename)
                    success = True
            # If not, check if available at server:
            if(success == False):
                check = urllib.urlopen(remoteurl)
                if(check.code == 200):
                    # File available, download! 
                    if(debug): printf('file %s available at GFS server at %s -> downloading'%(filename,datetime.datetime.now().time()))
                    urllib.urlretrieve(remoteurl,localfile)
                else:
                    # File not (yet) available, sleep a while and re-do the checks 
                    time.sleep(300)
        except:
            # Something weird happened. Sleep a bit, try again
            printf('weird exception on file %s: '%filename)
            printf(sys.exc_info()[0]) 
            time.sleep(30)

## Downloads the requested GFS data
# @param olga Pointer to object with OLGA settings
def downloadGFS(olga,islice):
    printf('Obtaining GFS data...')
    remotebase = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs'
    gfsdir     = '%s%04i%02i%02i/'%(olga.gfsDataRoot,olga.year,olga.month,olga.day) # where to save files
    dtGFS      = 3. # Time step of GFS input data [h]

    # Check if GFS data directory exists, if not create it
    if not os.path.exists(gfsdir):
        if(debug): printf('Making directory %s'%gfsdir)
        os.mkdir(gfsdir)

    # Calculate first and last hour to download
    t0 = olga.tstart + (islice+0) * olga.tslice
    t1 = olga.tstart + (islice+1) * olga.tslice
    nt = int((t1-t0) / dtGFS) + 1 

    # Loop over time steps
    downloads = []
    for t in range(nt):
        tact = int(t0 + t * dtGFS) # Forecast time
        loc = '%s.%04i%02i%02i%02i'%(remotebase,olga.year,olga.month,olga.day,olga.cycle) # Location at server
        if(res == '0.25'):
            fil = 'gfs.t%02iz.pgrb2.0p25.f%03i'%(olga.cycle,tact) # File at server
        elif(res == 'old'):
            fil = 'gfs.t%02iz.pgrb2f%02i'%(olga.cycle,tact) # File at server, old GFS naming
        else:
            sys.exit('Resolution %s invalid'%res)

        remoteurl = '%s/%s'%(loc,fil) # Path to file at server
        localfile = gfsdir + fil # Path to file on hard disk 

        # Start process:
        downloads.append(Process(target=downloadGFSFile, args=(localfile, remoteurl, fil,)))
        downloads[-1].start()

    # Join processes:
    for download in downloads:
        download.join() 

    printf('finished GFS at %s'%datetime.datetime.now().time())

## replace everything after '=' on line containing 'searchstring' with 'value' in 'filein' 
# @param filein Path to file
# @param searchstring String to search for
# @param value Value to replace
def replace(filein,searchstring,value):
    execute('sed -i -e "s/\(' +searchstring+ r'\).*/\1 = ' +value+ '/g" ' + filein)

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
    printf('updating namelists WPS and WRF...')

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
    printf('Running WPS at %s'%datetime.datetime.now().time())
    # Grr, we have to call the routines from the directory itself...
    os.chdir(olga.wpsRoot)

    if(olga.islice==0):
        # Cleanup stuff from previous day
        execute('rm GRIBFILE* >& /dev/null')
        execute('rm FILE* >& /dev/null')
        execute('rm met_em* >& /dev/null')

    # Each log will be saved in olgaLogs with yyyymmddhh added
    ss = olga.startstruct
    logappend = '%04i%02i%02i%02i_r%i'%(ss.year,ss.month,ss.day,ss.hour,olga.islice)

    # Make directory for the logs, if it doesn't exist
    if not os.path.exists(olga.olgaLogs):
        os.mkdir(olga.olgaLogs)

    # Run geogrid
    if(debug): printf('... WPS -> geogrid')
    execute('./geogrid.exe >& %sgeogrid.%s'%(olga.olgaLogs,logappend))
    #execute('./geogrid.exe')

    # Link the GFS data to the WPS directory
    gfsData = '%s%04i%02i%02i'%(olga.gfsDataRoot,olga.year,olga.month,olga.day)
    execute('./link_grib.csh '+gfsData+'/gfs*')

    # Temp for switching between old and new GFS
    # --------------------
    execute('rm Vtable >& /dev/null')
    if(res == '0.25'):
        execute('cp Vtable.GFS_0.25d Vtable')
    elif(res == 'old'):
        execute('cp Vtable.GFS_0.5d Vtable')
    # --------------------

    if(debug): printf('... WPS -> ungrib')

    # Run ungrib
    execute('./ungrib.exe >& %sungrib.%s'%(olga.olgaLogs,logappend))

    # Run metgrid
    if(debug): printf('... WPS -> metgrid')
    execute('./metgrid.exe >& %smetgrid.%s'%(olga.olgaLogs,logappend))

    printf('finished WPS at %s'%datetime.datetime.now().time())
    os.chdir(olga.domainRoot)

## Run the WRF steps
# @param olga Pointer to object with OLGA settings
def execWRF(olga):
    printf('Running WRF at %s'%datetime.datetime.now().time())
    # Grr, we have to call the routines from the directory itself...
    os.chdir(olga.wrfRoot)

    # Remove restart file and output from previous day(s)
    if(olga.islice==0):
        execute('rm wrfrst* >& /dev/null')
        execute('rm wrfout* >& /dev/null')

    # Each log will be saved in olgaLogs with yyyymmddhh added
    ss = olga.startstruct
    logappend = '%04i%02i%02i%02i_r%i'%(ss.year,ss.month,ss.day,ss.hour,olga.islice)

    # Make directory for the logs, if it doesn't exist
    if not os.path.exists(olga.olgaLogs):
        os.mkdir(olga.olgaLogs)
 
    # Link the met_em input files 
    execute('rm met_em*')
    execute('ln -s '+olga.wpsRoot+'met_em* .')

    #if(olga.ompThreads > 1):
    execute('export OMP_NUM_THREADS=%i'%(olga.ompThreads))
    execute('export | grep OMP')

    # Run real
    if(debug): printf('... WRF -> real.exe')
    if(olga.mpiTasks > 1):
        execute('mpirun -n %i ./real.exe >& %sreal.%s'%(olga.mpiTasks,olga.olgaLogs,logappend))
    else:
        execute('./real.exe >& %sreal.%s'%(olga.olgaLogs,logappend))

    # Run WRF as background process to allow other processes (download GFS, ..) to run at the same time..
    if(debug): printf('... WRF -> wrf.exe')
    if(olga.mpiTasks > 1):
        execute('mpirun -n %i ./wrf.exe >& %swrf.%s &'%(olga.mpiTasks,olga.olgaLogs,logappend))
    else:
        execute('./wrf.exe >& %swrf.%s &'%(olga.olgaLogs,logappend))

    # Load balancer Thunder
    #subprocess.call('sbatch run.slurm',shell=True,executable='/bin/bash')

    os.chdir(olga.domainRoot)

## Wait until the required restart file is available (i.e. WRF finished)
# @param olga Pointer to object with OLGA settings
def wait4WRF(olga):
    startTime = datetime.datetime.now()
    printf('Waiting for WRF to finish')

    es = olga.endstruct
    wrfrst = '%swrfrst_d01_%04i-%02i-%02i_%02i:%02i:00'%(olga.wrfRoot,es.year,es.month,es.day,es.hour,es.minute)

    # Check if 'wrfrst' is available, else sleep
    maxExecutionWRF = 10800  # Maximum time allowed for WRF, in case restart is never written
    while(True):
        elapsedTime = datetime.datetime.now()-startTime
        #printf(elapsedTime.seconds)
        if(elapsedTime.seconds > maxExecutionWRF):
            printf("WRF hit runtime limit of %i s!!!"%maxExecutionWRF)
            return 1 

        if(not os.path.isfile(wrfrst)):
            time.sleep(10)
        else: 
            printf('finished WRF at %s'%datetime.datetime.now().time())
            return 0

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
            execute('ncrcat -O %s %s%s'%(tmp,olga.wrfDataRoot,outname))

## Create the plots / maps / soundings
# @param olga Pointer to object with OLGA settings
def execPlots(olga):
    startTime = datetime.datetime.now()
    printf('Starting plots at %s'%startTime)

    # Spawn different processes for each mape type, saves quite some time..
    for dom in range(olga.ndom):
        if(olga.maps[dom]==True):
            printf('making maps for domain=%i'%(dom+1))
            #plot_driver(olga,dom,'maps')
            pmap = Process(target=plot_driver, args=(olga,dom,'maps',))
            pmap.start()

        if(olga.meteogr[dom]):
            printf('making time series for domain=%i'%(dom+1))
            #plot_driver(olga,dom,'time')
            pmgram = Process(target=plot_driver, args=(olga,dom,'time',))
            pmgram.start()

        if(olga.sounding[dom]==True):
            printf('making soundings for domain=%i'%(dom+1))
            #plot_driver(olga,dom,'sounding')
            psound = Process(target=plot_driver, args=(olga,dom,'sounding',))
            psound.start()

        if(olga.maps[dom]==True):
            pmap.join()
        if(olga.meteogr[dom]):
            pmgram.join()
        if(olga.sounding[dom]==True):
            psound.join()

    endTime = datetime.datetime.now()
    printf('finished plots at %s, execution took %s'%(endTime, endTime-startTime))

## Create the plots / maps / soundings
# @param olga Pointer to object with OLGA settings
def uploadPlots(olga):
    startTime = datetime.datetime.now()
    printf('Starting upload at %s'%startTime)
 
    # What's the best way to do this.....? For now hard coded
    local  = '%s%04i%02i%02i_t%02iz'%(olga.figRoot,olga.year,olga.month,olga.day,olga.cycle)
    remote = 'vanstratum-com@ssh.pcextreme.nl:~/domains/vanstratum.com/htdocs/olga/results/'
    execute('scp -r %s %s &> /dev/null'%(local, remote))

    endTime = datetime.datetime.now()
    printf('finished upload at %s'%endTime)
