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
#from netCDF4 import Dataset 
from copy import deepcopy
import sys
import datetime

from constants import *
from tools import *

def timestring2time(ts_in):
    return float("".join(ts_in[11:13])) + float("".join(ts_in[14:16]))/60.

def av(v1,v2):
    return (v1 + v2) * 0.5

"""
given day-of-year, time [utc], lat and lon, return potential (clear sky)
shortwave incoming radiation
"""
def get_swin(doy,time,lat,lon):
    doy    = float(doy)
    time   = float(time)
    lon    = -lon   #?! 
    sda    = 0.409 * np.cos(2. * np.pi * (doy - 173.) / 365.)
    sinlea = np.sin(2. * np.pi * lat / 360.) * np.sin(sda) - \
             np.cos(2. * np.pi * lat / 360.) * np.cos(sda) * \
             np.cos(2. * np.pi * (time*3600.) / 86400. - 2. * np.pi * lon / 360.)
    sinlea = max(sinlea, 0.0001)
    Tr     = (0.6 + 0.2 * sinlea)
    return 1368. * Tr * sinlea


"""
Given updraft height and velocity, calculate potential cross-country (constant height) velocity
Input can be nD arrays, or single values. see constants.py for definitions a,b,c,peff,etc.
"""
def VgemCrossCountry(z_upd,w_upd):
    vstf    = -(b-(b-(4.*-a*(-c-w_upd))**0.5))/(2.*-a) # speed to fly given updraft (MacCready) velocity [km h-1]
    wstf    = -a*vstf**2.+b*vstf-c # sink glider at vstf [m s-1]
    alpha   = -wstf / (w_upd - wstf) # fraction time spent circling [-]
    Vgem    = (1.-alpha)*vstf*peff # (1-alpha)*V = cross-country speed, correct for efficiency pilot [km h-1] 

    # Decrease potential speed for some arbitrary (to-do: tune) conditions
    if(np.size(z_upd)>1):
        Vgem[np.where(z_upd<800)] = Vgem[np.where(z_upd<800)] / 2.
        Vgem[np.where(z_upd<500)] = 0.
    else:
        Vgem = Vgem / 2. if z_upd<800 else Pv
        Vgem = 0.        if z_upd<500 else Pv

    return Vgem


"""
Read in data: all locations -> all time (mainly 2d fields)
"""
class readwrf_all:
    def __init__(self,file,domain):
        # Check if the netcdf4-python module is available:
        try:
            from netCDF4 import Dataset
            netcdf4py = True 
        except:
            netcdf4py = False 

        # Fallback option: NetCDF from Scientific.IO
        try:
            from Scientific.IO import NetCDF
            netcdfsci = True
        except:
            netcdfsci = False

        if(netcdf4py):
            print 'blaalt'
            postprocess_nc(self)
        elif(netcdfsci):
            read_nc_scientificIO(self,file,domain)
            postprocess_nc(self)
        else:
            sys.exit('No NETCDF module available')

def read_nc_scientificIO(obj,file,domain):
    print 'reading file %s'%file
    wrfin            = NetCDF.NetCDFFile(file,'r')
    nt               = len(wrfin.variables["HFX"][:,0,0])  

"""
Read in data: all locations -> all time (mainly 2d fields)
"""
class readwrf_all:
    def __init__(self,file,domain):
        # Check if the netcdf4-python module is available:
        try:
            from netCDF4 import Dataset
            wrfin = Dataset(file,'r')
            netcdf4py = True 
        except:
            netcdf4py = False 

        # Fallback option: NetCDF from Scientific.IO
        try:
            from Scientific.IO import NetCDF
            wrfin = NetCDF.NetCDFFile(file,'r')
            netcdfsci = True
        except:
            netcdfsci = False

        if(netcdf4py==False and netcdfsci==False):
            sys.exit('No NetCDF module available')

        print 'reading file %s'%file

        nt               = len(wrfin.variables["HFX"][:,0,0])  
        # In our case, {lat/lon/hgt} doesn't change in time since we don't have moving domains... 
        self.lat         = wrfin.variables["XLAT"][0,:,:] # latitude
        self.nlat        = np.size(self.lat[0,:])
        self.lon         = wrfin.variables["XLONG"][0,:,:] # longitude
        self.nlon        = np.size(self.lon[0,:])
        self.hgt         = wrfin.variables["HGT"][0,:,:] # terrain height 

        # Base state variables
        self.T00         = wrfin.variables["T00"][:] # reference temperature [K]
        self.P00         = wrfin.variables["P00"][:] # reference pressure [pa]
        self.ps          = wrfin.variables["PSFC"][:,:,:] # surface pressure [pa] 
        self.T2          = wrfin.variables["T2"][:,:,:] # 2m temperature [K]

        # read in for all domains:
        self.hfx         = wrfin.variables["HFX"][:,:,:] # sensible heat flux [W/m2]
        self.lh          = wrfin.variables["LH"][:,:,:] # latent heat flux [W/m2]
        self.rr_mp       = wrfin.variables["RAINNC"][:,:,:] # total microphysical rain [mm]
        self.rr_con      = wrfin.variables["RAINC"][:,:,:] # total convective rain [mm]
        self.U10         = wrfin.variables["U10"][:,:,:] # 10m u-wind [m/s]
        self.V10         = wrfin.variables["V10"][:,:,:] # 10m v-wind [m/s]
        self.U1000       = (wrfin.variables["U"][:,25,:,1:] + wrfin.variables["U"][:,25,:,:-1])/2.
        self.V1000       = (wrfin.variables["V"][:,25,1:,:] + wrfin.variables["V"][:,25,:-1,:])/2.
        self.ustar       = wrfin.variables["UST"][:,:,:] # surface friction velocity [m/s]
        self.slps        = wrfin.variables["PSFC"][:,:,:] \
                           / (1.-2.25577e-5 * self.hgt[:,:])**5.25588 # sea level pressure
        self.swd         = wrfin.variables["SWDOWN"][:,:,:] # shortwave incomming rad at sfc [W/m2]

        # REALLLLY ugly (and incorrect), but seems to work quite okay...:
        #   in theory: if one grid level cloud cover = 100%, total column should be 100%
        #   in WRF: this creates a 0% or 100% cloud cover switch. Averaging seems to do better.....
        #   to-do: weighted average? 
        self.cclow       = np.sum(wrfin.variables["CLDFRA"][:,:35,:,:],axis=1)   / 35. 
        self.ccmid       = np.sum(wrfin.variables["CLDFRA"][:,35:52,:,:],axis=1) / 17. 
        self.cchig       = np.sum(wrfin.variables["CLDFRA"][:,52:,:,:],axis=1)   / 11.
        self.ccsum       = self.cclow + self.ccmid + self.cchig
        self.ccsum[np.where(self.ccsum>1.)] = 1.

        # Get datetime in format "YYYY-MM-DD HH:MM:SS"
        self.times       = wrfin.variables["Times"][:,:] # date-time as individual characters
        self.time        = wrfin.variables["XTIME"][:] * 60. # time since start of simulation [s]
        self.dt          = timestring2time(self.times[1])-timestring2time(self.times[0])
        self.datetime    = [] # date-time as merged string
        self.date        = [] # date-time as merged string
        self.hour        = [] # hour of day
        for t in range(nt):
            self.datetime.append("".join(self.times[t,:10])+' '+"".join(self.times[t,11:19])) 
            self.date.append("".join(self.times[t,:10])) 
            self.hour.append(float(self.times[t,11] + self.times[t,12]))

        try: # Try if the TEMF (bl_pbl_physics=10) variables are available:
            self.zi      = wrfin.variables["HD_TEMF"][:,:,:] # dry thermal top TEMF
            self.zct     = wrfin.variables["HCT_TEMF"][:,:,:] # cloud top TEMF
            self.TEMF    = True
        except KeyError: # Use PBLH from other schemes, set cloud top to zero
            self.zi      = wrfin.variables["PBLH"][:,:,:] # boundary layer height [m]
            self.zct     = np.zeros_like(self.zi) + filval 
            self.TEMF    = False

        # Derived variables:
        rhos             = self.ps / (Rd * self.T2) # surface density [kg m-3]
        wthvs            = (self.hfx/(rhos*cp)) + 0.61*self.T2*(self.lh/(rhos*Lv)) # surface buoyancy flux [W m-2]
        wthvs[np.where(wthvs==0)] = eps # remove zero flux to prevent div/0
        self.L           = -(self.ustar**3. * tref) / (kappa * g * wthvs) # Obukhov length [m]
        wthvs[np.where(wthvs<0)] = 0. # remove negative flux for w* calculation 
        self.wstar       = (g * self.zi * wthvs / tref)**(1./3.) # convective velocity scale w* [m s-1]
        self.wglider     = deepcopy(self.wstar - supd) # w* minus sink glider [m s-1]
        self.wstar[np.where(self.wglider<0)] = 0. # set minimum updraft velocity to zero

        # Calculate the PFD per day 
        self.PFD = np.zeros((self.nlon,self.nlat))
        pV  = VgemCrossCountry(self.zi[:],self.wstar[:])
        # integrate to obtain cumulative flyable distance
        for t2 in range(1,pV[:,0,0].size):
            self.PFD[:,:] += pV[t2-1,:,:] * self.dt 

"""
Read in data from single location -> all time
"""
class readwrf_loc:
    def __init__(self,file,domain,glon,glat):
        # Check if the netcdf4-python module is available:
        try:
            from netCDF4 import Dataset
            wrfin = Dataset(file,'r')
            netcdf4py = True 
        except:
            netcdf4py = False 

        # Fallback option: NetCDF from Scientific.IO
        try:
            from Scientific.IO import NetCDF
            wrfin = NetCDF.NetCDFFile(file,'r')
            netcdfsci = True
        except:
            netcdfsci = False

        if(netcdf4py==False and netcdfsci==False):
            sys.exit('No NetCDF module available')

        n                = 3     # +/- how many grid points to average
        n1               = n+1 
        nt               = len(wrfin.variables["XTIME"][:]) # number of time steps in output
        self.lat         = wrfin.variables["XLAT"][0,:,:] # latitude 
        self.lon         = wrfin.variables["XLONG"][0,:,:] # longitude

        # Find gridpoint closest to glon,glat
        idx              = (((self.lat-glat)**2.+(self.lon-glon)**2.)**0.5).argmin()
        i1               = idx/float(len(self.lat[0,:]))
        jj               = int(np.floor(i1))
        ii               = int((i1-jj)*len(self.lat[0,:]))
        self.lat         = self.lat[jj,ii]
        self.lon         = self.lon[jj,ii]
        self.time        = wrfin.variables["XTIME"][:] * 60.
        self.T00         = wrfin.variables["T00"][:] # base state temperature [K]
        self.P00         = wrfin.variables["P00"][:] # base state pressure [pa]
        self.p           = wrfin.variables["P"][:,:,jj,ii] + wrfin.variables["PB"][:,:,jj,ii] # pressure [pa]
        self.phyd        = wrfin.variables["P_HYD"][:,:,jj,ii] # hydrostatic pressure [pa] (difference with pressure??)
        self.z           = (wrfin.variables["PH"][:,:,jj,ii] + wrfin.variables["PHB"][:,:,jj,ii]) / g # height [m]
        self.zf          = (self.z[:,1:]+self.z[:,:-1])/2. # height at full (mid gridpoint) levels [m]
        self.hgt         = wrfin.variables["HGT"][:,jj,ii] # terrain height 

        # Get datetime in format "YYYY-MM-DD HH:MM:SS"
        datetime         = wrfin.variables["Times"][:,:] # timedate array
        self.dt          = timestring2time(datetime[1])-timestring2time(datetime[0])
        self.datetime    = []
        self.date        = []
        self.hour        = [] # hour of day
        for t in range(nt):
            self.datetime.append("".join(datetime[t,:10])+' '+"".join(datetime[t,11:19])) 
            self.date.append("".join(datetime[t,:10])) 
            self.hour.append(float(datetime[t,11] + (datetime[t,12])))

        # For the next variables, average over area
        self.ps          = np.apply_over_axes(np.mean,wrfin.variables["PSFC"]  [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # surface pressure [pa]
        self.T2          = np.apply_over_axes(np.mean,wrfin.variables["T2"]    [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # 2m temperature [K]
        self.hfx         = np.apply_over_axes(np.mean,wrfin.variables["HFX"]   [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # sensible heat flux [W m-2]
        self.lh          = np.apply_over_axes(np.mean,wrfin.variables["LH"]    [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # latent heat flux [W m-2]
        self.th          = np.apply_over_axes(np.mean,wrfin.variables["T"]     [:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # potential temperature [K]
        self.qv          = np.apply_over_axes(np.mean,wrfin.variables["QVAPOR"][:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # vapor mixing ratio [kg kg-1]
        self.ql          = np.apply_over_axes(np.mean,wrfin.variables["QCLOUD"][:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # cloud liquid water [kg kg-1]
        self.cc          = np.apply_over_axes(np.mean,wrfin.variables["CLDFRA"][:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # cloud fraction [-]
        self.u           = np.apply_over_axes(np.mean,wrfin.variables["U"]     [:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # u-wind component [m s-1]
        self.v           = np.apply_over_axes(np.mean,wrfin.variables["V"]     [:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # v-wind component [m s-1]
        self.u10         = np.apply_over_axes(np.mean,wrfin.variables["U10"]   [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # 10m u-wind [m s-1] 
        self.v10         = np.apply_over_axes(np.mean,wrfin.variables["V10"]   [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # 10m v-wind [m s-1]
        self.q2          = np.apply_over_axes(np.mean,wrfin.variables["Q2"]    [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # 2m vapor mixing ratio [kg kg-1]

        self.ccl         = np.zeros((3,nt))
        self.ccl[0,:]    = np.apply_over_axes(np.mean,wrfin.variables["CLDFRA"][:,:35,  jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,0,0,0] 
        self.ccl[1,:]    = np.apply_over_axes(np.mean,wrfin.variables["CLDFRA"][:,35:52,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,0,0,0] 
        self.ccl[2,:]    = np.apply_over_axes(np.mean,wrfin.variables["CLDFRA"][:,52:,  jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,0,0,0] 

        try: # Try if the TEMF (bl_pbl_physics=10) variables are available:
            self.w       = np.apply_over_axes(np.mean,wrfin.variables["WUPD_TEMF"][:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # updraft velocity TEMF 
            self.zi      = np.apply_over_axes(np.mean,wrfin.variables["HD_TEMF"]  [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # dry thermal top TEMF
            self.ct      = np.apply_over_axes(np.mean,wrfin.variables["HCT_TEMF"] [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # cloud top TEMF
            self.lcl     = np.apply_over_axes(np.mean,wrfin.variables["LCL_TEMF"] [:,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # lifting condens. level TEMF
            self.thtemf  = np.apply_over_axes(np.mean,wrfin.variables["THUP_TEMF"][:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # potential temperature updrafts [K]
            self.qttemf  = np.apply_over_axes(np.mean,wrfin.variables["QTUP_TEMF"][:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # total water mixing rat. updrafts [kg kg-1]
            self.qltemf  = np.apply_over_axes(np.mean,wrfin.variables["QLUP_TEMF"][:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # liq. water mixing rat. updrafts [K]
            self.c3dtemf = np.apply_over_axes(np.mean,wrfin.variables["CF3D_TEMF"][:,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # cloud fraction updrafts [K]
            self.TEMF    = True
        except KeyError: # Other schemes don't have these variables, set to filval
            self.zi      = np.apply_over_axes(np.mean,wrfin.variables["PBLH"][:,jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0]
            self.w       = np.zeros_like(self.zi) + filval
            self.zi      = np.zeros_like(self.zi) + filval
            self.ct      = np.zeros_like(self.zi) + filval
            self.lcl     = np.zeros_like(self.zi) + filval
            self.thtemf  = np.zeros_like(self.zi) + filval
            self.qttemf  = np.zeros_like(self.zi) + filval
            self.qltemf  = np.zeros_like(self.zi) + filval
            self.c3dtemf = np.zeros_like(self.zi) + filval
            self.TEMF    = False

        # Derived variables:
        rhos             = self.ps / (Rd * self.T2) # surface density [kg m-3]
        wthvs            = (self.hfx/(rhos*cp)) + 0.61*self.T2*(self.lh/(rhos*Lv)) # surface buoyancy flux [W m-2]
        wthvs[np.where(wthvs<0)] = 0. # remove negative flux for w* calc 
        self.wstar       = (g * self.zi * wthvs/tref)**(1./3.) # convective velocity scale w*
        self.wstar[np.where(self.wstar<supd)] = 0. # convective velocity scale w* - sink glider
        self.pV          = VgemCrossCountry(self.zi,self.wstar) # potential cross-country velocity
        self.cPFD        = np.zeros_like(self.pV) # cumulative potential flight distance
        for t in range(1,nt):
            self.cPFD[t] = self.cPFD[t-1] + av(self.pV[t-1],self.pV[t]) * self.dt
            if(self.hour[t] < self.hour[t-1]): # new day: reset PFD
                self.cPFD[t] = 0.

        # Calculate temperature and dewpoint from potential temperature and total water mixing ratio
        self.qt  = self.qv + self.ql  # total water mixing ratio [kg kg-1]
        self.T   = np.zeros_like(self.th) # absolute temperature [K]
        self.Td  = np.zeros_like(self.th) # dew point temperature [K]
        self.Tu  = np.zeros_like(self.th) # absolute temperature updrafts [K]
        self.Tdu = np.zeros_like(self.th) # dew point temperature updrafts [K]

        for t in range(nt):
            self.th[t,:]  = self.th[t,:] + 300. # potential temp = base state (300) + perturbation (T)
            self.T[t,:]   = self.th[t,:] * (self.p[t,:] / 1.e5)**(287.05/1004.) # 1.e5 = reference pressure
            e             = ((self.p[t,:]) * self.qt[t,:]) / ((Rd/Rv) + self.qt[t,:]) # vapor pressure
            self.Td[t,:]  = ((1./273.) - (Rv/Lv) * np.log(e/611.))**-1.

            if(self.TEMF):
                self.Tu[t,:]  = self.thtemf[t,:] * (self.p[t,:] / 1.e5)**(287.05/1004.) # 1.e5 = reference pressure
                e             = ((self.p[t,:]) * self.qttemf[t,:]) / ((Rd/Rv) + self.qttemf[t,:]) # vapor pressure
                self.Tdu[t,:] = ((1./273.) - (Rv/Lv) * np.log(e/611.))**-1.
            else:
                self.Tu[t,:]  = -9999
                self.Tdu[t,:] = -9999
