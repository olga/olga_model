import datetime
import sys

# Import OLGA specific routines
from src.main import *

# Import file with domain settings from directory "domain_test" or other
from OLGA_NL.domainSettings import olgaSettings as settings_d1
olga = settings_d1()

olga.year      = 2015 
olga.month     = 02
olga.day       = 18

olga.maps      = ([False])
olga.sounding  = ([False])
olga.meteogr   = ([True ])

vars1          = (['pfd','wstar','zidry','cudepth','convection','swd','rain','wind10m','wind1000m'])         
olga.map_vars  = ([vars1])

print('--------------------------------')
print('Starting OLGA for %02i-%02i-%04i %02iz'%(olga.day, olga.month, olga.year, olga.cycle))
print('Start time: %s'%(datetime.datetime.now()))
print('--------------------------------')

# Loop over the time slices
startTime = datetime.datetime.now()
nslice = int(olga.ttotal/olga.tslice)
for islice in range(nslice):
    print('--------------------------------')
    print('Processing time slice %i of %i'%(islice+1,nslice))
    print('--------------------------------')

    olga.set_time(islice) # update time settings for the current time slice

    if(True):
        execPlots(olga)
    if(True):
        uploadPlots(olga) 

print('--------------------------------')
print('Execution time OLGA: %s'%(datetime.datetime.now()-startTime))
print('--------------------------------')
