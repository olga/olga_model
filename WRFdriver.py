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

from src.plotdriver import plot_driver

debug = True  

## Main settings for OLGA
# wrapped in an object to simplify passing the settings around
class olga_Settings:
    def __init__(self):
        # Local file system settings. APPEND EVERY PATH with '/'! 
        self.wpsRoot      = '/home/bart/WRFnl/WPSV3/' # Path to root of WPS
        self.wrfRoot      = '/home/bart/WRFnl/WRFV3/run_eu/' # Path to root of WRF
        self.olgaRoot     = '/home/bart/WRFnl/olga/' # Full path to OLGA scripts
        self.olgaLogs     = '/home/bart/WRFnl/olga/logs/' # Location to save logs
        self.figRoot      = '/home/scratch1/WRFnl/olga_results/' # Path to save OLGA figures
        self.gfsDataRoot  = '/home/scratch1/WRFnl/GFSdata/' # Path to store the GFS data
        self.wrfDataRoot  = '/home/scratch1/WRFnl/dataWRF/' # Path to store the WRF output

        # Computational settings. 
        self.mpi_tasks    = 2 # Number of MPI tasks
        self.omp_thr      = 2 # Number of OpenMP threads

        # Number of domains. This can be less than the size of the arrays below
        # in which case only the first ndom are used (e.g. for quick testing of outer domain)
        self.ndom         = 1  

        # Time settings
        self.ttotal       = 48 # Total time to simulate [h]
        self.tslice       = 24 # Split 'ttotal' in 'tslice' chunks [h]
        self.dt_output    = ([60,30]) # 'history_interval' from namelist, per domain, in minutes

        # Manually specify times (UTC) over which to make PFD's and time series. Only if both times are within
        # one 'tslice', maps are made. BvS: add better description :)
        # Same 'tanalysis' is used for all domains!
        self.tanalysis    = ([4,20])

        # Main map settings per domain
        self.maps         = ([True,True]) # Make soundings or not
        self.map_lat      = ([49.5, 51.3]) # Central latitude of map [deg]
        self.map_lon      = ([6.0,  6.7]) # Central longitude of map [deg]
        self.map_width    = ([1000000,590000]) # Domain plot width [m]
        self.map_height   = ([1000000,590000]) # Domain plot height [m]
        self.map_res      = (['l','i']) # Details of map (c=crude, l=low, i=interm, h=high)
        self.map_desc     = (['18x18km','6x6km']) 

        # Plot variables maps
        vars1 = (['pfd','swd','wstar','zidry','clouds','rr','wind10m','wind1000m'])         
        vars2 = (['pfd','wstar','zidry','cudepth'])         
        self.map_vars     = ([vars1,vars2]) # variables to plot per domain

        # Settings soundings (Detailed settings are in src/skewtlogp.py)
        self.sounding     = ([True,True]) # Make soundings or not
        sound_lat1        = ([self.map_lat[0]]) # Tmp arrays to populate locations, default=map center 
        sound_lon1        = ([self.map_lon[0]])  
        sound_lat2        = ([self.map_lat[1]]) 
        sound_lon2        = ([self.map_lon[1]]) 
        sound_name1       = (['center1'])
        sound_name2       = (['center2'])
        self.sound_name   = ([sound_name1,sound_name2]) # description of sounding location !! NO SPACES !! 
        self.sound_lat    = ([sound_lat1,sound_lat2]) # sounding latitudes 
        self.sound_lon    = ([sound_lon1,sound_lon2]) # sounding longitudes

        # Settings time series (meteograms)
        self.meteogr      = ([True,True]) # Make meteogram or not
        meteog_lat1       = ([self.map_lat[0]]) # Tmp arrays to populate locations, default=map center
        meteog_lon1       = ([self.map_lon[0]])  
        meteog_lat2       = ([self.map_lat[1]]) 
        meteog_lon2       = ([self.map_lon[1]])  
        meteog_name1      = (['center1'])
        meteog_name2      = (['center2'])
        self.meteogr_name = ([meteog_name1,meteog_name2]) # description of sounding location !! NO SPACES !! 
        self.meteogr_lat  = ([sound_lat1,sound_lat2]) # sounding latitudes 
        self.meteogr_lon  = ([sound_lon1,sound_lon2]) # sounding longitudes

        # -----------------------------------------------------
        # Don't change below, some checking of input
        if(self.ttotal%self.tslice != 0):
            sys.exit('ttotal should be integerer multiple of tslice')

    # Use time structs to do calculations on time  
    def set_time(self,islice):
        self.islice      = islice
        self.abs_start   = datetime.datetime.strptime('%i %i %i %i'%(self.day,self.month,self.year,self.tstart),"%d %m %Y %H")
        self.t0          = self.tstart + (self.islice+0) * self.tslice
        self.t1          = self.tstart + (self.islice+1) * self.tslice
        self.startstruct = self.abs_start + datetime.timedelta(hours=self.t0)
        self.endstruct   = self.abs_start + datetime.timedelta(hours=self.t1)

## Downloads the requested GFS data, or waits until available
# @param olga Pointer to object with OLGA settings
def download_GFS(olga,islice):
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
def update_Namelists(olga):
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
def run_WPS(olga):
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
    if(not os.path.isfile('Vtable')):
        subprocess.call('ln -s ungrib/Variable_Tables/Vtable.GFS Vtable',shell=True,executable='/bin/bash')

    # Run ungrib
    subprocess.call('./ungrib.exe >& %sungrib.%s'%(olga.olgaLogs,logappend),shell=True,executable='/bin/bash')

    # Run metgrid
    if(debug): print('... WPS -> metgrid')
    subprocess.call('./metgrid.exe >& %smetgrid.%s'%(olga.olgaLogs,logappend),shell=True,executable='/bin/bash')

    print('finished WPS at %s'%datetime.datetime.now().time())
    os.chdir(olga.olgaRoot)

## Run the WRF steps
# @param olga Pointer to object with OLGA settings
def run_WRF(olga):
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

    if(olga.omp_thr > 1):
        subprocess.call('export OMP_NUM_THREADS=%i'%(olga.omp_thr),shell=True,executable='/bin/bash')

    # Run real
    if(debug): print('... WRF -> real.exe')
    if(olga.mpi_tasks > 1):
        subprocess.call('mpirun -n %i ./real.exe >& %sreal.%s'%(olga.mpi_tasks,olga.olgaLogs,logappend),shell=True,executable='/bin/bash')
    else:
        subprocess.call('./real.exe >& %sreal.%s'%(olga.olgaLogs,logappend),shell=True,executable='/bin/bash')

    # Run WRF as background process to allow other processes (download GFS, ..) to run at the same time..
    if(debug): print('... WRF -> wrf.exe')
    if(olga.mpi_tasks > 1):
        subprocess.call('mpirun -n %i ./wrf.exe >& %swrf.%s &'%(olga.mpi_tasks,olga.olgaLogs,logappend),shell=True,executable='/bin/bash')
    else:
        subprocess.call('./wrf.exe >& %swrf.%s &'%(olga.olgaLogs,logappend),shell=True,executable='/bin/bash')

    os.chdir(olga.olgaRoot)

## Wait until the required restart file is available (i.e. WRF finished)
# @param olga Pointer to object with OLGA settings
def wait_for_WRF(olga):
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
def move_WRF_output(olga):
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
def make_plots(olga):
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

        pmap.join()
        pmgram.join()
        psound.join()

    print('finished plots at %s'%datetime.datetime.now().time())

## Main OLGA function, called when calling WRFdriver.py
# @param mode ....
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
    year   = 2014 #time.strftime('%Y')
    month  = 06   #time.strftime('%m')
    day    = 26   #time.strftime('%d')
    tstart = 00   # start time of simulation 
    cycle  = 0    # which GFS cycle? {0,6,12,18}

    # Create object with OLGA settings, and add time settings
    olga_settings         = olga_Settings()
    olga_settings.year    = int(year)
    olga_settings.month   = int(month)
    olga_settings.day     = int(day)
    olga_settings.tstart  = int(tstart)
    olga_settings.cycle   = int(cycle)
    olga_settings.basestr = "%s%s%s_t%02iz"%(year,month,day,cycle) # base string for output

    # Loop over the time slices
    startTime = datetime.datetime.now()
    nslice = int(olga_settings.ttotal/olga_settings.tslice)

    print('--------------------------------')
    print('Starting OLGA: %s'%(datetime.datetime.now()))
    print('--------------------------------')

    # Loop over the different time slices
    for islice in range(nslice):
        print('--------------------------------')
        print('Processing time slice %i of %i'%(islice+1,nslice))
        print('--------------------------------')

        olga_settings.set_time(islice) # update time settings for the current time slice

        if(mode == 'all'):
            download_GFS(olga_settings,islice) # download GFS data
            update_Namelists(olga_settings) # update WRf & WPS namelists
            run_WPS(olga_settings) # run the WPS routines
            run_WRF(olga_settings) # run the WRF routines
            wait_for_WRF(olga_settings) # Wait until the restart file is available
        if(mode == 'all' or mode == 'post'):
            move_WRF_output(olga_settings)
            make_plots(olga_settings)

    print('--------------------------------')
    print('Execution time OLGA: %s'%(datetime.datetime.now()-startTime))
    print('--------------------------------')


