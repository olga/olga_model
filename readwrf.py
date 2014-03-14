from constants import *
import numpy as np
from netCDF4 import Dataset 
from copy import deepcopy
import sys
from pylab import *

# ---------------------------------------------
# Calculate PFD from MacCready theory
#----------------------------------------------
class getpfd:
  def __init__(self,file):

    wrfin = Dataset(file,'r')
    nt    = len(wrfin.variables["HFX"][:,0,0])
    nlat  = np.size(wrfin.variables["XLAT"][0,0,:] )
    nlon  = np.size(wrfin.variables["XLONG"][0,0,:])
    hfx   = wrfin.variables["HFX"][:,:,:]
    lh    = wrfin.variables["LH"][:,:,:]
    ps    = wrfin.variables["PSFC"][:,:,:] 
    T2    = wrfin.variables["T2"][:,:,:] 

    # Find the time step in WRF output
    times = wrfin.variables["Times"][:]
    t0 = float("".join(times[0][11:13])) + float("".join(times[0][14:16]))/60.
    t1 = float("".join(times[1][11:13])) + float("".join(times[1][14:16]))/60.
    dth = t1-t0
    print "dt in PFD calc = %f"%dth

    # Get (dry) updraft height from TEMF or other PBL scheme
    try:
      pblh           = wrfin.variables["HD_TEMF"][:,:,:]
    except KeyError:
      try:
        pblh         = wrfin.variables["PBLH"][:,:,:]
      except KeyError:
        print "found no pblh..."

    # Calculate surface buoyancy flux
    rhos = ps / (Rd * T2)                           # surface density
    wthvs = (hfx/(rhos*cp)) + 0.61*T2*(lh/(rhos*Lv)) # surface buoy. flux
    wthvs[np.where(wthvs<0)] = 0.                    # remove negative flux for w* calc 
    wstar = (g * pblh * wthvs/tref)**(1./3.)         # w* (convective velocity scale)

    # Set wstar to zero where wstar < sink glider
    wstar[np.where(wstar<supd)] = 0.

    # Set wstar to zero if updraft height less than 500m
    wstar[np.where(pblh<500.)] = 0.

    vstf    = -(b-(b-(4.*-a*(-c-wstar))**0.5))/(2.*-a)
    wstf    = -a*vstf**2.+b*vstf-c
    alpha   = -wstf / (wstar - wstf)
    self.pV = (1.-alpha)*vstf*peff  

    # Decrease potential speed if updraft height less then 800m
    self.pV[np.where(pblh<800)] = self.pV[np.where(pblh<800)] / 2.

    self.pfd         = np.zeros((nt,nlon,nlat))
    for t in range(1,nt):
      self.pfd[t,:,:] = self.pfd[t-1,:,:] + (self.pV[t-1,:,:] * dth)

# ---------------------------------------------
# Read in data: all locations -> all time (mainly 2d fields)
#----------------------------------------------
class readwrf_all:
  def __init__(self,file,domain):
    print 'reading file %s'%file

    wrfin            = Dataset(file,'r')
    nt               = len(wrfin.variables["HFX"][:,0,0])  
 
    # In our case, {lat/lon/hgt} doesn't change in time since we don't have moving domains... 
    self.lat         = wrfin.variables["XLAT"][0,:,:]  ; self.nlat = np.size(self.lat[0,:])
    self.lon         = wrfin.variables["XLONG"][0,:,:] ; self.nlon = np.size(self.lon[0,:])
    self.hgt         = wrfin.variables["HGT"][0,:,:]          # terrain height 

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
    self.ustar       = wrfin.variables["UST"][:,:,:]          # 10m v-wind [m/s]
    self.slps        = wrfin.variables["PSFC"][:,:,:] / (1.-2.25577e-5 * self.hgt[:,:])**5.25588

    # REALLLLY ugly (and incorrect), but seems to work quit okay...:
    #   in theory: if one grid level cloud cover = 100%, total column should be 100%
    #   in WRF: this creates a 0% or 100% cloud cover switch. Averaging seems to do better.....
    #   to-do: weighted average? 
    self.cclow       = np.sum(wrfin.variables["CLDFRA"][:,:35,:,:],axis=1)   / 35. 
    self.ccmid       = np.sum(wrfin.variables["CLDFRA"][:,35:52,:,:],axis=1) / 17. 
    self.cchig       = np.sum(wrfin.variables["CLDFRA"][:,52:,:,:],axis=1)   / 11.
    self.ccsum       = self.cclow + self.ccmid + self.cchig
    self.ccsum[np.where(self.ccsum>1.)] = 1.

    # Synthethic albedo
    # Get very crude estimate of density profile for LWP integration
    self.thvref      = (wrfin.variables["T"][0,:,0,0] + self.T00[0]) * (1.-(1.-Rv/Rd) * wrfin.variables["QVAPOR"][0,:,0,0])
    self.pref        = wrfin.variables["P"][0,:,0,0] + wrfin.variables["PB"][0,:,0,0]
    self.dnref       = self.pref / (Rd * self.thvref)
    self.z           = (wrfin.variables["PH"][0,:,0,0] + wrfin.variables["PHB"][0,:,0,0]) / g
    self.zf          = (self.z[1:]+self.z[:-1])/2. 
    self.dz          = self.z[1:] - self.z[:-1]

    self.LWP         = np.sum(wrfin.variables["QCLOUD"][:,:,:,:] * self.dnref[None,:,None,None] * self.dz[None,:,None,None],axis=1)
    self.IWP         = np.sum(wrfin.variables["QICE"][:,:,:,:] * self.dnref[None,:,None,None] * self.dz[None,:,None,None],axis=1)
    self.CWP         = self.LWP + self.IWP
    cnm = 400.
    tau = 0.19 * self.CWP**(5./6.)*cnm**(1./3.) 
    self.calb        = tau / (6.8 + tau) 

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
    self.L = -(self.ustar**3. * tref) / (kappa * g * wthvs)    # !Obukhov length
    wthvs[np.where(wthvs<0)] = 0.                             # remove negative flux for w* calc 
    self.wstar = (g * self.zi * wthvs/tref)**(1./3.)
    self.wstar[np.where(self.wstar<supd)] = 0.                # convective velocity scale w* - sink glider


# ---------------------------------------------
# Read in data from single location -> all time
#----------------------------------------------
class readwrf_loc:
  def __init__(self,file,domain,glon,glat):
    import sys
   
    # +/- how many grid points to average
    n  = 3   
    n1 = n+1 

    wrfin            = Dataset(file,'r')
    nt               = len(wrfin.variables["XTIME"][:])
    self.lat         = wrfin.variables["XLAT"][0,:,:]  
    self.lon         = wrfin.variables["XLONG"][0,:,:] 

    # Find gridpoint nearest to glon,glat
    idx              = (((self.lat-glat)**2.+(self.lon-glon)**2.)**0.5).argmin()
    i1               = idx/float(len(self.lat[0,:]))
    jj             = int(np.floor(i1))
    ii             = int((i1-jj)*len(self.lat[0,:]))
    self.lat         = self.lat[jj,ii]
    self.lon         = self.lon[jj,ii]

    print 'Reading 1D at lon=%.3f, lat=%.3f'%(self.lon,self.lat)

    self.time        = wrfin.variables["XTIME"][:] * 60.
    self.T00         = wrfin.variables["T00"][:]
    self.P00         = wrfin.variables["P00"][:]
    self.p           = wrfin.variables["P"][:,:,jj,ii] + wrfin.variables["PB"][:,:,jj,ii]
    self.phyd        = wrfin.variables["P_HYD"][:,:,jj,ii]
    self.z           = (wrfin.variables["PH"][:,:,jj,ii] + wrfin.variables["PHB"][:,:,jj,ii]) / g
    self.zf          = (self.z[:,1:]+self.z[:,:-1])/2. 

    # Get date-time and merge chars to string
    datetime         = wrfin.variables["Times"][:,:]          # timedate array
    self.datetime    = []
    # Get datetime in format "YYYY-MM-DD HH:MM:SS"
    for t in range(nt):
      self.datetime.append("".join(datetime[t,:10])+' '+"".join(datetime[t,11:19])) 

    # For the next variables, average over area
    self.ps          = np.apply_over_axes(np.mean,wrfin.variables["PSFC"]  [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0]
    self.T2          = np.apply_over_axes(np.mean,wrfin.variables["T2"]    [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] 
    self.hfx         = np.apply_over_axes(np.mean,wrfin.variables["HFX"]   [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] 
    self.lh          = np.apply_over_axes(np.mean,wrfin.variables["LH"]    [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] 
    self.th          = np.apply_over_axes(np.mean,wrfin.variables["T"]     [:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0]
    self.qv          = np.apply_over_axes(np.mean,wrfin.variables["QVAPOR"][:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0]
    self.ql          = np.apply_over_axes(np.mean,wrfin.variables["QCLOUD"][:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0]
    self.cc          = np.apply_over_axes(np.mean,wrfin.variables["CLDFRA"][:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0]
    self.u           = np.apply_over_axes(np.mean,wrfin.variables["U"]     [:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0]
    self.v           = np.apply_over_axes(np.mean,wrfin.variables["V"]     [:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0]
    self.u10         = np.apply_over_axes(np.mean,wrfin.variables["U10"]   [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0]
    self.v10         = np.apply_over_axes(np.mean,wrfin.variables["V10"]   [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0]
    self.T2          = np.apply_over_axes(np.mean,wrfin.variables["T2"]    [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0]
    self.q2          = np.apply_over_axes(np.mean,wrfin.variables["Q2"]    [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0]

    self.ccl         = np.zeros((3,nt))
    self.ccl[0,:]    = np.apply_over_axes(np.mean,wrfin.variables["CLDFRA"][:,:35,  jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,0,0,0] 
    self.ccl[1,:]    = np.apply_over_axes(np.mean,wrfin.variables["CLDFRA"][:,35:52,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,0,0,0] 
    self.ccl[2,:]    = np.apply_over_axes(np.mean,wrfin.variables["CLDFRA"][:,52:,  jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,0,0,0] 

    if(domain==2): 
      self.w         = np.apply_over_axes(np.mean,wrfin.variables["WUPD_TEMF"][:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0]  # updraft velocity 
      self.zi        = np.apply_over_axes(np.mean,wrfin.variables["HD_TEMF"]  [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0]  # dry thermal top TEMF
      self.ct        = np.apply_over_axes(np.mean,wrfin.variables["HCT_TEMF"] [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0]  # cloud top TEMF

      # TEST: average updraft velocity from TEMF over ABL depth 
      self.wav = np.zeros(nt)
      for t in range(nt):
        if(self.zi[t] > 500):  
          k300    = (np.abs(self.zf[t,:]-300)).argmin()
          kzi     = (np.abs(self.zf[t,:]-self.zi[t])).argmin()
          self.wav[t]  = np.average(self.w[t,k300:kzi])

    elif(domain==1):
      self.zi        = np.apply_over_axes(np.mean,wrfin.variables["PBLH"][:,jj-n:jj+n1,ii-n:ii+n1][1,2])[:,0,0]

    rhos  = self.ps / (Rd * self.T2)
    wthvs = (self.hfx/(rhos*cp)) + 0.61*self.T2*(self.lh/(rhos*Lv))
    wthvs[np.where(wthvs<0)] = 0.                             # remove negative flux for w* calc 
    self.wstar = (g * self.zi * wthvs/tref)**(1./3.)
    self.wstar[np.where(self.wstar<supd)] = 0.               # convective velocity scale w* - sink glider

    self.qt          = self.qv + self.ql
    self.t           = np.zeros_like(self.th)
    self.pf          = np.zeros_like(self.th)
    self.T           = np.zeros_like(self.th)
    self.Td          = np.zeros_like(self.th)

    for t in range(nt):
      self.th[t,:]   = self.th[t,:] + 300. # Total potential temp = base state (300) + perturbation (T)
      self.T[t,:]    = self.th[t,:] * (self.p[t,:] / 1.e5)**(287.05/1004.) # 1.e5 = reference pressure
      e              = ((self.p[t,:]) * self.qt[t,:]) / ((Rd/Rv) + self.qt[t,:])
      self.Td[t,:]   = ((1./273.) - (Rv/Lv) * np.log(e/611.))**-1.

