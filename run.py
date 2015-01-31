import datetime

# Import OLGA specific routines
from src.main import *

# Import file with domain settings from directory "domain_test" or other
#from domain_test.domainSettings import olgaSettings as settings_d1
from OLGA_NL.domainSettings import olgaSettings as settings_d1

mode = 'post'

olga = settings_d1()

n = -1
olga.year    = 2014 #int(time.strftime('%Y'))
olga.month   = 6    #int(time.strftime('%m'))
olga.day     = 26   #int(time.strftime('%d'))
olga.tstart  = 00   # start time of simulation
olga.cycle   = 0    # which GFS cycle? {0,6,12,18}

olga.basestr = "%s%s%s_t%02iz"%(olga.year, olga.month, olga.day, olga.cycle) # base string for output

print('--------------------------------')
print('Starting OLGA: %s'%(datetime.datetime.now()))
print('--------------------------------')

# Loop over the time slices
startTime = datetime.datetime.now()
nslice = int(olga.ttotal/olga.tslice)
for islice in range(nslice):
    print('--------------------------------')
    print('Processing time slice %i of %i'%(islice+1,nslice))
    print('--------------------------------')

    olga.set_time(islice) # update time settings for the current time slice

    if(mode == 'all'):
        downloadGFS(olga,islice) # download GFS data
        updateNamelists(olga) # update WRf & WPS namelists
        execWPS(olga) # run the WPS routines
        execWRF(olga) # run the WRF routines
        wait4WRF(olga) # Wait until the restart file is available
        moveWRFOutput(olga)
    if(mode == 'all' or mode == 'post'):
        execPlots(olga)

print('--------------------------------')
print('Execution time OLGA: %s'%(datetime.datetime.now()-startTime))
print('--------------------------------')
