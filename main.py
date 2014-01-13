from plot import *

# 1. get command line arguments (Y,M,D,dom)
# --------------------------------------------------------
if(len(sys.argv) != 5):
  sys.exit('provide input, e.g. "2013 06 02 2" (year month day domain)')
else:
  year   = int(sys.argv[1])
  month  = int(sys.argv[2])
  day    = int(sys.argv[3])
  dom    = int(sys.argv[4])
  smooth = False

# 2. Check if wrfout available
# --------------------------------------------------------
date   = str(year)+str(month).zfill(2)+str(day).zfill(2)
wrfout = '../dataWRF/'+date+'/wrfout_d0'+str(dom)+'_'+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_00:00:00' 
if(not os.path.isfile(wrfout)):
  sys.exit('file %s not available'%wrfout) 

# --------------------------------------------------------
# 3. Create image directory if not available:
if not os.path.exists('figures/'+date):
  os.makedirs('figures/'+date)

# 4. Settings:
# --------------------------------------------------------
t0 = 12                                         # first time map
t1 = 14                                        # last time map
if(dom==1):
  variables = (['pfd','wstar','zidry','clouds','rr']) 
elif(dom==2):
  #variables = (['pfd','wstar','zidry','cudepth']) 
  variables = (['wstar']) 

# 5. Setup basemap only once:
# --------------------------------------------------------
if(dom==1):
  m = Basemap(width=1800000,height=1800000,
              rsphere=(6378137.00,6356752.3142),\
              resolution='l',area_thresh=10.,projection='lcc',\
              lat_1=51.4,lat_2=51.4,lat_0=51.4,lon_0=5.5)

elif(dom==2):
  m = Basemap(width=590000,height=590000,
              rsphere=(6378137.00,6356752.3142),\
              resolution='i',area_thresh=10.,projection='lcc',\
              lat_1=51.3,lat_2=51.3,lat_0=51.3,lon_0=6.7)

# --------------------------------------------------------
# 6. Create maps
create_maps(wrfout,dom,date,t0,t1,variables,m,filter=True)

  
