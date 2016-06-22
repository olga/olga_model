import datetime
import math

# Import OLGA specific routines
from src.main import *

# Import file with domain settings from directory "domain_test" or other
#from domain_test.domainSettings import olgaSettings as settings_d1
from OLGA_NL.domainSettings import olgaSettings as settings_d1

def myfloor(x, base=6):
	return int(base * math.floor(float(x)/base))

mode = 'all'

olga = settings_d1()

olga.year    = int(time.strftime('%Y'))
olga.month   = int(time.strftime('%m'))
olga.day     = int(time.strftime('%d'))
#olga.cycle   = myfloor(int(time.strftime('%H')))

print('--------------------------------')
print('Starting OLGA for %02i-%02i-%04i %02iz'%(olga.day, olga.month, olga.year, olga.cycle))
print('Start time: %s'%(datetime.datetime.now()))
print('--------------------------------')

execute('rm -rf %s/*' % olga.figRoot)
execute('rm -rf %s/*' % olga.wrfDataRoot)
execute('rm -rf %s/*' % olga.gfsDataRoot)
execute('rm -rf %s/*' % olga.olgaLogs)



# Loop over the time slices
startTime = datetime.datetime.now()
nslice = int(olga.ttotal/olga.tslice)
for islice in range(0, nslice):
    print('--------------------------------')
    print('Processing time slice %i of %i'%(islice+1,nslice))
    print('--------------------------------')

    olga.set_time(islice) # update time settings for the current time slice

    nerr = 0
    if(mode == 'all'):
        downloadGFS(olga,islice) # download GFS data
        updateNamelists(olga) # update WRf & WPS namelists
        execWPS(olga) # run the WPS routines
        nerr = execWRF(olga) # run the WRF routines

#        nerr = wait4WRF(olga) # Wait until the restart file is available
        if(nerr > 0):
            print("Something went wrong with WRF.... Stopping after this time slice")

    if(mode == 'all' or mode == 'post'):
        moveWRFOutput(olga) # Merge NetCDF output, and move to output dir
        execPlots(olga) # Make all maps, time series, etc.
    if(mode == 'all'):
        uploadPlots(olga) # Upload to server

    if(nerr > 0):
        break

local = olga.gfsDataRoot
remote = "wrf@cumulus-humilis.nl:~/"
execute("scp -l 8192 -r %s %s" % (local, remote))

print('--------------------------------')
print('Execution time OLGA: %s'%(datetime.datetime.now()-startTime))
print('--------------------------------')
