from constants import *
import numpy as np
from netCDF4 import Dataset 
from copy import deepcopy

# ---------------------------------------------
# Calculate PFD from MacCready theory
#----------------------------------------------
class getpfd:
  def __init__(self,file):
    wrfin            = Dataset(file,'r')
    nt               = len(wrfin.variables["HFX"][:,0,0])
    nlat             = np.size(wrfin.variables["XLAT"][0,:] )
    nlon             = np.size(wrfin.variables["XLONG"][0,:])
    hfx              = wrfin.variables["HFX"][:,:,:]
    lh               = wrfin.variables["LH"][:,:,:]
    ps               = wrfin.variables["PSFC"][:,:,:] 
    T2               = wrfin.variables["T2"][:,:,:] 

    # Get (dry) updraft height from TEMF or other PBL scheme
    try:
      pblh           = wrfin.variables["HD_TEMF"][:,:,:]
    except KeyError:
      try:
        pblh         = wrfin.variables["PBLH"][:,:,:]
      except KeyError:
        print "found no pblh..."

    # Calculate surface buoyancy flux
    rhos  = ps / (Rd * T2)                           # surface density
    wthvs = (hfx/(rhos*cp)) + 0.61*T2*(lh/(rhos*Lv)) # surface buoy. flux
    wthvs[np.where(wthvs<0)] = 0.                    # remove negative flux for w* calc 
    wstar = (g * pblh * wthvs/tref)**(1./3.)         # w* (convective velocity scale)
    wstar[np.where(wstar<supd)] = 0.                 # remove sink glider

    vstf             = -(b-(b-(4.*-a*(-c-wstar))**0.5))/(2.*-a)
    wstf             = -a*vstf**2.+b*vstf-c
    alpha            = -wstf / (wstar - wstf)
    self.pV          = (1.-alpha)*vstf*peff  

    self.pfd         = np.zeros((nt,nlon,nlat))
    for t in range(1,nt):
      self.pfd[t,:,:] = self.pfd[t-1,:,:] + self.pV[t-1,:,:]



# ---------------------------------------------
# Read in data: all locations -> all time
#----------------------------------------------
class readwrf_all:
  def __init__(self,file,domain):
    print 'reading file %s'%file

    wrfin            = Dataset(file,'r')
    nt               = len(wrfin.variables["HFX"][:,0,0])  
 
    # In our case, {lat/lon/hgt} doesn't change in time since we don't have moving domains... 
    self.lat         = wrfin.variables["XLAT"][:,:]  ; self.nlat = np.size(self.lat[0,:])
    self.lon         = wrfin.variables["XLONG"][:,:] ; self.nlon = np.size(self.lon[0,:])
    self.hgt         = wrfin.variables["HGT"][:,:]          # terrain height 

    # Base state variables
    self.T00         = wrfin.variables["T00"][:]
    self.P00         = wrfin.variables["P00"][:]
    self.ps          = wrfin.variables["PSFC"][:,:,:] 
    self.T2          = wrfin.variables["T2"][:,:,:]             # 2m temperature [K]

    # read in for all domains:
    self.hfx         = wrfin.variables["HFX"][:,:,:]          # sensible heat flux [W/m2]
    self.lh          = wrfin.variables["LH"][:,:,:]           # latent heat flux [W/m2]
    self.rr_mp       = wrfin.variables["RAINNC"][:,:,:]       # total microphysical rain [mm]
    self.rr_con      = wrfin.variables["RAINC"][:,:,:]        # total convective rain [mm]
    self.U10         = wrfin.variables["U10"][:,:,:]          # 10m u-wind [m/s]
    self.V10         = wrfin.variables["V10"][:,:,:]          # 10m v-wind [m/s]
    self.slps        = wrfin.variables["PSFC"][:,:,:] / (1.-2.25577e-5 * self.hgt[:,:])**5.25588

    # REALLLLY ugly (and incorrect), but seems to work quit okay...:
    #   in theory: if one grid level cloud cover = 100%, total column should be 100%
    #   in WRF: this creates a 0% or 100% cloud cover switch. Averaging seems to do better.....
    #   to-do: weighted average? 
    self.cclow       = wrfin.variables["CLDFRA"][:,0,:,:] 
    self.ccmid       = wrfin.variables["CLDFRA"][:,1,:,:] 
    self.cchig       = wrfin.variables["CLDFRA"][:,2,:,:] 
    self.ccsum       = self.cclow + self.ccmid + self.cchig
    self.ccsum[np.where(self.ccsum>1.)] = 1.

    # Get date-time and merge chars to string
    datetime         = wrfin.variables["Times"][:,:]          # timedate array
    self.datetime    = []

    # Get datetime in format "YYYY-MM-DD HH:MM:SS"
    for t in range(nt):
      self.datetime.append("".join(datetime[t,:10])+' '+"".join(datetime[t,11:19])) 

    # specific for domain 1 (large):
    if(domain==1):
      self.zi        = wrfin.variables["PBLH"][:,:,:]         # boundary layer height [m]

    # specific for domain 2 (small):
    elif(domain==2):
      self.zi        = wrfin.variables["HD_TEMF"][:,:,:]      # dry thermal top TEMF
      self.zct       = wrfin.variables["HCT_TEMF"][:,:,:]     # cloud top TEMF

    # Derived variables:
    rhos  = self.ps / (Rd * self.T2)
    wthvs = (self.hfx/(rhos*cp)) + 0.61*self.T2*(self.lh/(rhos*Lv))
    wthvs[np.where(wthvs<0)] = 0.                             # remove negative flux for w* calc 
    self.wstar = (g * self.zi[:,:] * wthvs/tref)**(1./3.)
    self.wstar[np.where(self.wstar<supd)] = 0.               # convective velocity scale w* - sink glider









#class getpfd:
#  def __init__(self,file):
#    wrfin            = Dataset(file,'r')
#    nt               = len(wrfin.variables["XTIME"][:])
#    nlat             = np.size(wrfin.variables["XLAT"][0,0,:] )
#    nlon             = np.size(wrfin.variables["XLONG"][0,0,:])
#
#    hfx              = wrfin.variables["HFX"][:,:,:]
#    hfx[np.where(hfx < 0.)] = 0.
#    pblh             = wrfin.variables["PBLH"][:,:,:]
#    wstar            = (9.81 * pblh * (hfx / (1.2*1004.)) / 300.)**0.333
#    wstar[np.where(wstar < supd)] = 0.
#
#    vstf             = -(b-(b-(4.*-a*(-c-wstar))**0.5))/(2.*-a)
#    wstf             = -a*vstf**2.+b*vstf-c
#    alpha            = -wstf / (wstar - wstf)
#    self.pV          = (1.-alpha)*vstf*peff  
#
#    self.pfd         = np.zeros((nt,nlon,nlat))
#    for t in range(1,nt):
#      self.pfd[t,:,:] = self.pfd[t-1,:,:] + self.pV[t-1,:,:]




# ---------------------------------------------
# Read in data from single location -> all time
#----------------------------------------------
#class readwrf_loc:
#  def __init__(self,file,domain,glon,glat):
#    import sys
#
#    wrfin            = Dataset(file,'r')
#    nt               = len(wrfin.variables["XTIME"][:])
#    self.lat         = wrfin.variables["XLAT"][0,:,:]  
#    self.lon         = wrfin.variables["XLONG"][0,:,:] 
#
#    # Find gridpoint nearest to glon,glat
#    idx              = (((self.lat-glat)**2.+(self.lon-glon)**2.)**0.5).argmin()
#    i1               = idx/float(len(self.lat[0,:]))
#    glat             = int(np.floor(i1))
#    glon             = int((i1-glat)*len(self.lat[0,:]))
#    self.lat         = self.lat[glat,glon]
#    self.lon         = self.lon[glat,glon]
#
#    print 'Reading 1D at lon=%.3f, lat=%.3f'%(self.lon,self.lat)
#
#    self.time        = wrfin.variables["XTIME"][:] * 60.
#    self.T00         = wrfin.variables["T00"][:]
#    self.P00         = wrfin.variables["P00"][:]
#    self.p           = wrfin.variables["P"][:,:,glat,glon] + wrfin.variables["PB"][:,:,glat,glon]
#    self.z           = (wrfin.variables["PH"][:,:,glat,glon] + wrfin.variables["PHB"][:,:,glat,glon]) / g
#    self.zf          = (self.z[:,1:]+self.z[:,:-1])/2. 
#
#    self.ps          = wrfin.variables["PSFC"][:,glat,glon] 
#    self.T2          = wrfin.variables["T2"][:,glat,glon] 
#    self.hfx         = wrfin.variables["HFX"][:,glat,glon] 
#    self.lh          = wrfin.variables["LH"][:,glat,glon] 
#
#    self.th          = wrfin.variables["T"][:,:,glat,glon]
#    self.qv          = wrfin.variables["QVAPOR"][:,:,glat,glon]
#    self.ql          = wrfin.variables["QCLOUD"][:,:,glat,glon]
#    self.cc          = wrfin.variables["CLDFRA"][:,:,glat,glon]
#
#    if(domain==2): 
#      self.w         = wrfin.variables["WUPD_TEMF"][:,:,glat,glon] # updraft velocity 
#      self.zi        = wrfin.variables["HD_TEMF"][:,glat,glon]     # dry thermal top TEMF
#      self.ct        = wrfin.variables["HCT_TEMF"][:,glat,glon]    # cloud top TEMF
#
#      # TEST: average updraft velocity from TEMF over ABL depth 
#      self.wav = np.zeros(nt)
#      for t in range(nt):
#        if(self.zi[t] > 500):  
#          k300    = (np.abs(self.zf[t,:]-300)).argmin()
#          kzi     = (np.abs(self.zf[t,:]-self.zi[t])).argmin()
#          self.wav[t]  = np.average(self.w[t,k300:kzi])
#
#    elif(domain==1):
#      self.zi        = wrfin.variables["PBLH"][:,glat,glon]
#
#    rhos  = self.ps / (Rd * self.T2)
#    wthvs = (self.hfx/(rhos*cp)) + 0.61*self.T2*(self.lh/(rhos*Lv))
#    wthvs[np.where(wthvs<0)] = 0.                             # remove negative flux for w* calc 
#    self.wstar = (g * self.zi * wthvs/tref)**(1./3.)
#    self.wstar[np.where(self.wstar<supd)] = 0.               # convective velocity scale w* - sink glider
#
#    #self.qt          = self.qv + self.ql
#    #self.t           = np.zeros_like(self.th)
#    #self.pf          = np.zeros_like(self.th)
#    #self.T           = np.zeros_like(self.th)
#    #self.Td          = np.zeros_like(self.th)
#    #self.Thu         = np.zeros_like(self.th)
#    #self.Tup         = np.zeros_like(self.th)
#
#    # CHECK CALCULATIONS!!!!!!!!!
#    #for t in range(nt):
#    #  self.th[t,:]   = self.th[t,:] + 300.
#    #  self.T[t,:]    = self.th[t,:] * (self.p[t,:] / 1.e5)**(287.05/1004.)
#    #  self.Td[t,:]   = (5.42e3 / np.log((0.622 * 2.53e11) / (self.qt[t,:] * self.p[t,:])))
#    #  if(temf):
#    #    self.Thu[t,:]  = self.thlTEMF[t,:] + (2.45e6/1004.)*self.ql[t,:]
#    #    self.Tup[t,:]  = self.Thu[t,:] * (self.p[t,:] / 1.e5)**(287.05/1004.)
#
##d = readwrf_loc('../dataWRF/20080512/wrfout_d02_2008-05-12_00:00:00',2,6.93146,52.2913)
#
#
#
#d = readwrf_all('../dataWRF/20140308/wrfout_d01_2014-03-08_00:00:00_2d.nc',domain=1)


# ---------------------------------------------
# Read in data: all locations -> all time
#----------------------------------------------
#class readwrf_all:
#  def __init__(self,file,domain):
#    print 'reading file %s'%file
#
#    wrfin            = Dataset(file,'r')
#    nt               = len(wrfin.variables["XTIME"][:])
#   
#    # In our case, {lat/lon/hgt} doesn't change in time since we don't have moving domains... 
#    self.lat         = wrfin.variables["XLAT"][0,:,:]  ; self.nlat = np.size(self.lat[0,:])
#    self.lon         = wrfin.variables["XLONG"][0,:,:] ; self.nlon = np.size(self.lon[0,:])
#    self.hgt         = wrfin.variables["HGT"][0,:,:]          # terrain height 
#
#    # Base state variables
#    self.T00         = wrfin.variables["T00"][:]
#    self.P00         = wrfin.variables["P00"][:]
#    self.ps          = wrfin.variables["PSFC"][:,:,:] 
#    self.T2          = wrfin.variables["T2"][:,:,:]             # 2m temperature [K]
#
#    # read in for all domains:
#    self.hfx         = wrfin.variables["HFX"][:,:,:]          # sensible heat flux [W/m2]
#    self.lh          = wrfin.variables["LH"][:,:,:]           # latent heat flux [W/m2]
#    self.rr_mp       = wrfin.variables["RAINNC"][:,:,:]       # total microphysical rain [mm]
#    self.rr_con      = wrfin.variables["RAINC"][:,:,:]        # total convective rain [mm]
#
#    # REALLLLY ugly (and incorrect), but seems to work quit okay...:
#    #   in theory: if one grid level cloud cover = 100%, total column should be 100%
#    #   in WRF: this creates a 0% or 100% cloud cover switch. Averaging seems to do better.....
#    #   to-do: weighted average? 
#    self.cclow       = np.sum(wrfin.variables["CLDFRA"][:,:35,:,:],axis=1)   / 35. 
#    self.ccmid       = np.sum(wrfin.variables["CLDFRA"][:,35:52,:,:],axis=1) / 17. 
#    self.cchig       = np.sum(wrfin.variables["CLDFRA"][:,52:,:,:],axis=1)   / 11.
#    self.ccsum       = self.cclow + self.ccmid + self.cchig
#    self.ccsum[np.where(self.ccsum>1.)] = 1.
#
#    # Get date-time and merge chars to string
#    datetime         = wrfin.variables["Times"][:,:]          # timedate array
#    self.datetime    = []
#
#    # Get datetime in format "YYYY-MM-DD HH:MM:SS"
#    for t in range(nt-2):
#      print "BvS hacked time!!"
#      self.datetime.append("".join(datetime[t,:10])+' '+"".join(datetime[t,11:19])) 
#
#    # specific for domain 1 (large):
#    if(domain==1):
#      self.zi        = wrfin.variables["PBLH"][:,:,:]         # boundary layer height [m]
#      self.U10       = wrfin.variables["U10"][:,:,:]          # 10m u-wind [m/s]
#      self.V10       = wrfin.variables["V10"][:,:,:]          # 10m v-wind [m/s]
#      self.slps      = wrfin.variables["PSFC"][:,:,:] / (1.-2.25577e-5 * self.hgt[:,:])**5.25588
#
#    # specific for domain 2 (small):
#    elif(domain==2):
#      self.zi        = wrfin.variables["HD_TEMF"][:,:,:]      # dry thermal top TEMF
#      self.zct       = wrfin.variables["HCT_TEMF"][:,:,:]     # cloud top TEMF
#
#    # Derived variables:
#    rhos  = self.ps / (Rd * self.T2)
#    wthvs = (self.hfx/(rhos*cp)) + 0.61*self.T2*(self.lh/(rhos*Lv))
#    wthvs[np.where(wthvs<0)] = 0.                             # remove negative flux for w* calc 
#    self.wstar = (g * self.zi[:,:] * wthvs/tref)**(1./3.)
#    self.wstar[np.where(self.wstar<-supd)] = 0.               # convective velocity scale w* - sink glider
 

