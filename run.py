import datetime
import math
import importlib

# Import OLGA specific routines
from src.main import *

#Usage: domain year month day cycle


# Import file with domain settings from directory "domain_test" or other
#from domain_test.domainSettings import olgaSettings as settings_d1



def myfloor(x, base=6):
	return int(base * math.floor(float(x)/base))

if len(sys.argv) < 2:
    print >>sys.stderr, "Usage: %s domain [mode] [year] [month] [day]" % sys.argv[0]
    sys.exit(-1)


domain = sys.argv[1]
domainSettings = importlib.import_module("%s.domainSettings" % domain)
settings_d1 = domainSettings.olgaSettings
olga = settings_d1()

mode = 'all'
olga.year    = int(time.strftime('%Y'))
olga.month   = int(time.strftime('%m'))
olga.day     = int(time.strftime('%d'))
#olga.cycle   = myfloor(int(time.strftime('%H')))

if len(sys.argv) > 2:
    mode = sys.argv[2]
if len(sys.argv) > 3:
    olga.year = int(sys.argv[3])
if len(sys.argv) > 4:
    olga.month = int(sys.argv[4])
if len(sys.argv) > 5:
    olga.day = int(sys.argv[5])

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
    if mode == 'all' or mode == 'download':
        downloadGFS(olga,islice) # download GFS data
   
    if mode == 'all':
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
remote = "wrf@baardman.net:~/"
execute("scp -l 8192 -r %s %s" % (local, remote))

print('--------------------------------')
print('Execution time OLGA: %s'%(datetime.datetime.now()-startTime))
print('--------------------------------')
