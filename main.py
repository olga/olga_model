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
#create_maps(wrfout,dom,date,t0,t1,variables,m,filter=True)


#locations  = ['de Bilt','Essen','Beauvecchain','Bergen','Idar Oberstein']
#longitudes = [ 5.18    ,  6.96 ,  4.77        ,  9.93  ,  7.33          ]
#latitudes  = [52.10    , 51.40 , 50.75        , 52.81  ,  49.70         ]
#
#close('all')
#
#def modplot(ax):
#  from matplotlib.ticker import AutoMinorLocator
#  #minorLocator   = AutoMinorLocator()
#  #ax.yaxis.set_minor_locator(minorLocator)
#  minorLocator   = AutoMinorLocator()
#  ax.xaxis.set_minor_locator(minorLocator)
#
#  ax.spines['right'].set_visible(False)
#  ax.get_yaxis().tick_left()
#  ax.spines['top'].set_visible(False)
#  ax.get_xaxis().tick_bottom()
#
#from colormaps import *
#
#d = readwrf_loc(wrfout,dom,6.96,51.40)
#figure()
#
#ax=subplot(511)
#contourf(d.time/3600.,d.zf[0,:],transpose(d.cc),cmap=wup)
##colorbar()
#xlim(0,24)
#modplot(ax)
#
#ax=subplot(512)
#plot(d.time/3600.,d.T2[:]-273.)
#xlim(0,24)
#modplot(ax)
#
#ax=subplot(513)
#bar(d.time/3600.,d.zi,width=0.2,label='zi TEMF')
#xlim(0,24)
#modplot(ax)
#
#ax=subplot(514)
#bar(d.time/3600.-0.1,d.wstar,color='b',width=0.2,label='wstar')
#bar(d.time/3600.+0.1,d.wav,color='g',width=0.2,label='w TEMF')
#legend()
#xlim(0,24)
#modplot(ax)
#
#ax=subplot(515)
#plot(d.time/3600.,d.hfx,label='sens. hf')
#plot(d.time/3600.,d.lh,label='lat. hf')
#legend()
#xlim(0,24)
#modplot(ax)




