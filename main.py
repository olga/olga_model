from plot import *

# --------------------------------------------------------
# 1. get command line arguments (Y,M,D,dom)
# --------------------------------------------------------
if(len(sys.argv) != 6):
  sys.exit('provide input, e.g. "2013 06 02 2 maps" (year month day domain type)')
else:
  year   = int(sys.argv[1])
  month  = int(sys.argv[2])
  day    = int(sys.argv[3])
  dom    = int(sys.argv[4])
  type   = sys.argv[5]
  dt     = 3
  smooth = False

# --------------------------------------------------------
# 2. Check if wrfout available
# --------------------------------------------------------
date   = str(year)+str(month).zfill(2)+str(day).zfill(2)
wrfout = '../dataWRF/'+date+'/wrfout_d0'+str(dom)+'_'+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_00:00:00' 
if(not os.path.isfile(wrfout)):
  sys.exit('file %s not available'%wrfout) 

# --------------------------------------------------------
# 3. Create image directory if not available
# --------------------------------------------------------
if not os.path.exists('figures/'+date):
  os.makedirs('figures/'+date)

# --------------------------------------------------------
# 4. Create maps
# --------------------------------------------------------
if(type=="maps"):
  # 1. Setup basemap only once and (deep)copy later
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

  t0 = 0     # first time map
  t1 = 24    # last time map
  if(dom==1):
    variables = (['pfd','clouds','rr','slpwind']) 
    #variables = (['zidry']) 
  elif(dom==2):
    #variables = (['pfd','wstar','zidry','cudepth']) 
    #variables = (['pfd','wstar']) 
    #variables = (['clouds']) 
    variables = (['zol']) 

  create_maps(wrfout,dom,date,t0,t1,dt,variables,m,filter=True)

# --------------------------------------------------------
# 5. Create meteograms
# --------------------------------------------------------
if(type=="tser"):
  names = []
  lats = []
  lons = []

  #names.append("Sterksel")
  #lats.append(51.35) 
  #lons.append(05.61)

  names.append("Aachen")
  lats.append(50.78) 
  lons.append(06.09)

  if(False):
    names.append("Deelen")
    lats.append(52.06) 
    lons.append(05.88)

    names.append("Antwerp")
    lats.append(51.22) 
    lons.append(4.41)

    names.append("Groningen")
    lats.append(53.22) 
    lons.append(6.55)

    names.append("Bremen")
    lats.append(53.08) 
    lons.append(8.82)

    names.append("Kassel")
    lats.append(51.32) 
    lons.append(9.465)

    names.append("Frankfurt")
    lats.append(50.12) 
    lons.append(8.69)

    names.append("Florennes")
    lats.append(50.25) 
    lons.append(4.61)

    names.append("Trier")
    lats.append(49.75) 
    lons.append(6.62)

    names.append("Meschede")
    lats.append(51.35) 
    lons.append(8.277)

  create_tser(wrfout,dom,date,names,lons,lats)


