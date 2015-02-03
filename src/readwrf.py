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
from copy import deepcopy
import sys
import datetime

from constants import *
from tools import *

# TMP BVS
from pylab import *

## Convert the output from WRF ("Times" variable) to a time in hours
# @param ts_in Single line of WRF's "Times" variable
def timestring2time(ts_in):
    return float("".join(ts_in[11:13])) + float("".join(ts_in[14:16]))/60.

## Some date-time conversions from WRF's "Times" variable to more useful merged format
# @param ncobj Pointer to the read NetCDF object
def wrf_time_conversion(ncobj):
    ncobj.dt          = timestring2time(ncobj.times[1])-timestring2time(ncobj.times[0]) # time step of model outut [h]
    ncobj.datetime    = [] # date-time as merged string
    ncobj.date        = [] # date as merged string
    ncobj.year        = [] # year as merged string
    ncobj.month       = [] # month as merged string
    ncobj.day         = [] # day as merged string
    ncobj.hour        = [] # hour of day as merged string
    ncobj.doy         = [] # day of the year
    for t in range(ncobj.nt):
        ncobj.datetime.append("".join(ncobj.times[t,:10])+' '+"".join(ncobj.times[t,11:19])) 
        ncobj.date.append("".join(ncobj.times[t,:10])) 
        ncobj.year.append(float("".join(ncobj.times[t,:4])))
        ncobj.month.append(float("".join(ncobj.times[t,5:7])))
        ncobj.day.append(float("".join(ncobj.times[t,8:10])))
        ncobj.hour.append(float(ncobj.times[t,11] + ncobj.times[t,12]))
        ncobj.doy.append(datetime.datetime.strptime('%i %i %i %i'%(ncobj.day[-1],\
                        ncobj.month[-1],ncobj.year[-1],ncobj.hour[-1]),"%d %m %Y %H").timetuple().tm_yday)

## Find the indices of the times in "tanalysis" in the current NetCDF file
# @param olga Pointer to the OLGA settings object
# @param ncobj Pointer to the read NetCDF object
def find_analysis_times(olga,ncobj):
    ncobj.t0_ana = [] ; t0 = -1
    ncobj.t1_ana = []
    for t in range(ncobj.nt-1):
        if(ncobj.hour[t] <= olga.tanalysis[0] and ncobj.hour[t+1] > olga.tanalysis[0]):
            t0 = t
        if(ncobj.hour[t] <= olga.tanalysis[1] and ncobj.hour[t+1] > olga.tanalysis[1] and t0 != -1):
            ncobj.t0_ana.append(t0)
            ncobj.t1_ana.append(t)
            t0 = -1

## Calculate average over two values
def av(v1,v2):
    return (v1 + v2) * 0.5

"""
given day-of-year, time [utc], lat and lon, return potential (clear sky)
shortwave incoming radiation
"""
def swin(doy, time, lat, lon, returnElev=False):
    doy    = float(doy)
    time   = float(time)
    lon    = -lon # ?! 
    sda    = 0.409 * np.cos(2. * np.pi * (doy - 173.) / 365.)
    sinlea = np.sin(2. * np.pi * lat / 360.) * np.sin(sda) - \
             np.cos(2. * np.pi * lat / 360.) * np.cos(sda) * \
             np.cos(2. * np.pi * (time*3600.) / 86400. - 2. * np.pi * lon / 360.)
    if(np.size(lat)>1):
        sinlea[np.where(sinlea <= 0)] = 1e-9
    else:
        sinlea = max(sinlea, 1e-9)
    Tr     = (0.6 + 0.2 * sinlea)
    swin   = 1368. * Tr * sinlea

    if(returnElev == True):
        return (swin, sinlea)
    else:
        return swin

"""
Given updraft height and velocity, calculate potential cross-country (constant height) velocity
Input can be nD arrays, or single values.
"""
def VgemCrossCountry(z_upd, w_upd, a, b, c, pdfEff):
    vstf    = -(b-(b-(4.*-a*(-c-w_upd))**0.5))/(2.*-a) # speed to fly given updraft (MacCready) velocity [km h-1]
    wstf    = -a*vstf**2.+b*vstf-c # sink glider at vstf [m s-1]
    alpha   = -wstf / (w_upd - wstf) # fraction time spent circling [-]
    Vgem    = (1.-alpha) * vstf * pdfEff # (1-alpha)*V = cross-country speed, correct for efficiency pilot [km h-1] 

    # Decrease potential speed for some arbitrary (to-do: tune) conditions
    if(np.size(z_upd)>1):
        Vgem[z_upd < 800] *= 0.5
        Vgem[z_upd < 500]  = 0.
    else:
        Vgem = Vgem * 0.5 if z_upd < 800 else Pv
        Vgem = 0.         if z_upd < 500 else Pv

    return Vgem

## Function to try different NetCDF python modules
# @param file Path to a NetCDF file
def get_nc_obj(file):
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
        sys.exit('No NetCDF module available!')
    else:
        return wrfin
    
"""
Read in data: all locations -> all time (mainly 2d fields)
Note: after reading, data is stored as [time, lat, lon]
"""
class readwrf_all:
    def __init__(self,olga,file_in,t0,t1):
        #print 'reading file %s'%file
        t1              += 1
        wrfin            = get_nc_obj(file_in)
        ntglob           = wrfin.variables["XTIME"][:].size # number of time steps in output
        nt               = wrfin.variables["XTIME"][t0:t1].size ; self.nt = nt 

        # In our case, {lat/lon/hgt} doesn't change in time since we don't have moving domains... 
        self.lat         = wrfin.variables["XLAT"][0,:,:] # latitude
        self.nlat        = np.size(self.lat[:,0])
        self.lon         = wrfin.variables["XLONG"][0,:,:] # longitude
        self.nlon        = np.size(self.lon[0,:])
        self.hgt         = wrfin.variables["HGT"][0,:,:] # terrain height 

        # Base state variables
        self.T00         = wrfin.variables["T00"][t0:t1] # reference temperature [K]
        self.P00         = wrfin.variables["P00"][t0:t1] # reference pressure [pa]
        self.ps          = wrfin.variables["PSFC"][t0:t1,:,:] # surface pressure [pa] 
        self.T2          = wrfin.variables["T2"][t0:t1,:,:] # 2m temperature [K]

        # read in for all domains:
        self.hfx         = wrfin.variables["HFX"][t0:t1,:,:] # sensible heat flux [W/m2]
        self.lh          = wrfin.variables["LH"][t0:t1,:,:] # latent heat flux [W/m2]
        self.rr_mp       = wrfin.variables["RAINNC"][t0:t1,:,:] # total microphysical rain [mm]
        self.rr_con      = wrfin.variables["RAINC"][t0:t1,:,:] # total convective rain [mm]
        self.U10         = wrfin.variables["U10"][t0:t1,:,:] # 10m u-wind [m/s]
        self.V10         = wrfin.variables["V10"][t0:t1,:,:] # 10m v-wind [m/s]
        self.U1000       = (wrfin.variables["U"][t0:t1,25,:,1:] + wrfin.variables["U"][t0:t1,25,:,:-1])/2.
        self.V1000       = (wrfin.variables["V"][t0:t1,25,1:,:] + wrfin.variables["V"][t0:t1,25,:-1,:])/2.
        self.ustar       = wrfin.variables["UST"][t0:t1,:,:] # surface friction velocity [m/s]
        self.slps        = wrfin.variables["PSFC"][t0:t1,:,:] \
                           / (1.-2.25577e-5 * self.hgt[:,:])**5.25588 # sea level pressure
        self.swd         = wrfin.variables["SWDNB"][t0:t1,:,:] # shortwave incomming radiation at surface [W/m2]
        self.swdc        = wrfin.variables["SWDNBC"][t0:t1,:,:] # clear sky shortwave incomming radiation at surface [W/m2]
        self.swdc[self.swdc==0] = 1e-9
 
        # Try if the TEMF (bl_pbl_physics=10) variables are available:
        try: 
            self.zi      = wrfin.variables["HD_TEMF"][t0:t1,:,:] # dry thermal top TEMF
            self.zct     = wrfin.variables["HCT_TEMF"][t0:t1,:,:] # cloud top TEMF
            self.TEMF    = True
        except KeyError: # Use PBLH from other schemes, set cloud top to zero
            self.zi      = wrfin.variables["PBLH"][t0:t1,:,:] # boundary layer height [m]
            self.zct     = np.zeros_like(self.zi) + filval 
            self.TEMF    = False

        # REALLLLY ugly (and incorrect), but seems to work quite okay...:
        #   in theory: if one grid level cloud cover = 100%, total column should be 100%
        #   in WRF: this creates a 0% or 100% cloud cover switch. Averaging seems to do better.....
        #   to-do: weighted average? 
        self.cclow       = np.sum(wrfin.variables["CLDFRA"][t0:t1,:35,:,:],axis=1)   / 35. 
        self.ccmid       = np.sum(wrfin.variables["CLDFRA"][t0:t1,35:52,:,:],axis=1) / 17. 
        self.cchig       = np.sum(wrfin.variables["CLDFRA"][t0:t1,52:,:,:],axis=1)   / 11.
        self.ccsum       = self.cclow + self.ccmid + self.cchig
        self.ccsum[np.where(self.ccsum>1.)] = 1.

        # Read in the different time variables, and do some conversions in "wrf_time_conversion"
        self.times       = wrfin.variables["Times"][t0:t1,:] # date-time as individual characters
        self.time        = wrfin.variables["XTIME"][t0:t1] * 60. # time since start of simulation [s]
        wrf_time_conversion(self)
        find_analysis_times(olga,self)

        # Potential & fractional incoming shortwave radiation. Masked at -1 for sunrise/set
        #self.swd_t       = np.zeros((self.nt, self.nlat, self.nlon))
        #self.swd_frac    = np.zeros((self.nt, self.nlat, self.nlon))
        #self.sun_elev    = np.zeros((self.nt, self.nlat, self.nlon))
        #for t in range(self.nt):
        #    self.swd_t[t,:,:], self.sun_elev[t,:,:] = swin(self.doy[t], self.hour[t], self.lat, self.lon, returnElev=True)
        #self.swd_frac    = self.swd / self.swd_t
        #self.swd_frac[self.sun_elev < 0.05] = -1
        self.swdf         = self.swd / self.swdc # Fraction of incoming shortwave radiation
        self.swdf[self.swdc < 1] = -1 # Mask for night 

        # Updraft velocity
        rhos             = self.ps / (Rd * self.T2) # surface density [kg m-3]
        wthvs            = (self.hfx/(rhos*cp)) + 0.61*self.T2*(self.lh/(rhos*Lv)) # surface buoyancy flux [W m-2]
        wthvs[np.where(wthvs==0)] = eps # remove zero flux to prevent div/0
        self.L           = -(self.ustar**3. * tref) / (kappa * g * wthvs) # Obukhov length [m]
        wthvs[np.where(wthvs<0)] = 0. # remove negative flux for w* calculation 
        self.wstar       = (g * self.zi * wthvs / tref)**(1./3.) # convective velocity scale w* [m s-1]
        self.wglider     = deepcopy(self.wstar) - olga.sinkGlider # w* minus sink glider [m s-1]
        self.wglider[self.wglider<0] = 0. # Limit updraft velocity glider to zero

        # Calculate the PFD, only if there is a full day.
        if(np.size(self.t0_ana) > 0):
            self.date_PFD = [] # List to store the datetime string of PFD calculation
            self.PFD = np.zeros((np.size(self.t0_ana), np.size(olga.pfdNames), self.nlat, self.nlon)) # Empty array for PFD

            # Loop over different days 
            for i in range(np.size(self.t0_ana)):
                t0 = self.t0_ana[i]
                t1 = self.t1_ana[i]+1
                self.date_PFD.append(self.date[t0])

                # Loop over different glider types
                for j in range(np.size(olga.pfdNames)):
                    # Calculate the instantaneous achievable cross-country velocity given updraft velocity. Some corrections are applied.
                    pV  = VgemCrossCountry(self.zi[t0:t1], self.wglider[t0:t1], olga.pfdA[j], olga.pfdB[j], olga.pfdC[j], olga.pfdEff[j])

                    # integrate to obtain cumulative flyable distance. 
                    # At each time, the average distance over the past ouput period is added
                    for t2 in range(1, pV.shape[0]):
                        self.PFD[i,j,:,:] += 0.5*(pV[t2-1,:,:]+pV[t2,:,:]) * self.dt 
        else:
            self.PFD = False

"""
Read in data from single location -> all time
"""
class readwrf_loc:
    def __init__(self,olga,file_in,glon,glat,t0,t1):
        t1 += 1
        n = 3 ; n1 = n+1
        wrfin            = get_nc_obj(file_in)
        ntglob           = wrfin.variables["XTIME"][:].size # number of time steps in output
        nt               = wrfin.variables["XTIME"][t0:t1].size ; self.nt = nt 
        self.lat         = wrfin.variables["XLAT"][0,:,:] # latitude 
        self.lon         = wrfin.variables["XLONG"][0,:,:] # longitude

        # Find gridpoint closest to glon,glat
        idx              = (((self.lat-glat)**2.+(self.lon-glon)**2.)**0.5).argmin()
        i1               = idx/float(len(self.lat[0,:]))
        jj               = int(np.floor(i1))
        ii               = int((i1-jj)*len(self.lat[0,:]))
        self.lat         = self.lat[jj,ii]
        self.lon         = self.lon[jj,ii]

        self.time        = wrfin.variables["XTIME"][t0:t1] * 60.
        self.T00         = wrfin.variables["T00"][t0:t1] # base state temperature [K]
        self.P00         = wrfin.variables["P00"][t0:t1] # base state pressure [pa]
        self.p           = wrfin.variables["P"][t0:t1,:,jj,ii] + wrfin.variables["PB"][t0:t1,:,jj,ii] # pressure [pa]
        self.phyd        = wrfin.variables["P_HYD"][t0:t1,:,jj,ii] # hydrostatic pressure [pa] (difference with pressure??)
        self.z           = (wrfin.variables["PH"][t0:t1,:,jj,ii] + wrfin.variables["PHB"][t0:t1,:,jj,ii]) / g # height [m]
        self.zf          = (self.z[t0:t1,1:]+self.z[t0:t1,:-1])/2. # height at full (mid gridpoint) levels [m]
        self.hgt         = wrfin.variables["HGT"][t0:t1,jj,ii] # terrain height 

        # Read in the different time variables, and do some conversions in "wrf_time_conversion"
        self.times       = wrfin.variables["Times"][t0:t1,:] # date-time as individual characters
        self.time        = wrfin.variables["XTIME"][t0:t1] * 60. # time since start of simulation [s]
        wrf_time_conversion(self)
        find_analysis_times(olga,self)

        # For the next variables, average over area
        self.ps          = np.apply_over_axes(np.mean,wrfin.variables["PSFC"]  [t0:t1,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # surface pressure [pa]
        self.T2          = np.apply_over_axes(np.mean,wrfin.variables["T2"]    [t0:t1,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # 2m temperature [K]
        self.hfx         = np.apply_over_axes(np.mean,wrfin.variables["HFX"]   [t0:t1,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # sensible heat flux [W m-2]
        self.lh          = np.apply_over_axes(np.mean,wrfin.variables["LH"]    [t0:t1,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # latent heat flux [W m-2]
        self.th          = np.apply_over_axes(np.mean,wrfin.variables["T"]     [t0:t1,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # potential temperature [K]
        self.qv          = np.apply_over_axes(np.mean,wrfin.variables["QVAPOR"][t0:t1,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # vapor mixing ratio [kg kg-1]
        self.ql          = np.apply_over_axes(np.mean,wrfin.variables["QCLOUD"][t0:t1,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # cloud liquid water [kg kg-1]
        self.cc          = np.apply_over_axes(np.mean,wrfin.variables["CLDFRA"][t0:t1,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # cloud fraction [-]
        self.u           = np.apply_over_axes(np.mean,wrfin.variables["U"]     [t0:t1,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # u-wind component [m s-1]
        self.v           = np.apply_over_axes(np.mean,wrfin.variables["V"]     [t0:t1,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # v-wind component [m s-1]
        self.u10         = np.apply_over_axes(np.mean,wrfin.variables["U10"]   [t0:t1,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # 10m u-wind [m s-1] 
        self.v10         = np.apply_over_axes(np.mean,wrfin.variables["V10"]   [t0:t1,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # 10m v-wind [m s-1]
        self.q2          = np.apply_over_axes(np.mean,wrfin.variables["Q2"]    [t0:t1,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # 2m vapor mixing ratio [kg kg-1]
        self.swd         = np.apply_over_axes(np.mean,wrfin.variables["SWDOWN"][t0:t1,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # shortwave incomming radiation at surface [W/m2]
        self.rr_mp       = np.apply_over_axes(np.mean,wrfin.variables["RAINNC"][t0:t1,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # total microphysical rain [mm]
        self.rr_con      = np.apply_over_axes(np.mean,wrfin.variables["RAINC"] [t0:t1,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # total convective rain [mm]

        # Rain over last period
        self.drr_mp      = np.zeros(self.rr_mp.size)
        self.drr_con     = np.zeros(self.rr_con.size)
        self.drr_mp[1:]  = self.rr_mp[1:] - self.rr_mp[:-1]
        self.drr_con[1:] = self.rr_con[1:] - self.rr_con[:-1]

        self.ccl         = np.zeros((3,nt))
        self.ccl[0,:]    = np.apply_over_axes(np.mean,wrfin.variables["CLDFRA"][t0:t1,:35,  jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,0,0,0] 
        self.ccl[1,:]    = np.apply_over_axes(np.mean,wrfin.variables["CLDFRA"][t0:t1,35:52,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,0,0,0] 
        self.ccl[2,:]    = np.apply_over_axes(np.mean,wrfin.variables["CLDFRA"][t0:t1,52:,  jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,0,0,0] 


        try: # Try if the TEMF (bl_pbl_physics=10) variables are available:
            self.w       = np.apply_over_axes(np.mean,wrfin.variables["WUPD_TEMF"][t0:t1,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # updraft velocity TEMF 
            self.zi      = np.apply_over_axes(np.mean,wrfin.variables["HD_TEMF"]  [t0:t1,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # dry thermal top TEMF
            self.ct      = np.apply_over_axes(np.mean,wrfin.variables["HCT_TEMF"] [t0:t1,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # cloud top TEMF
            self.lcl     = np.apply_over_axes(np.mean,wrfin.variables["LCL_TEMF"] [t0:t1,  jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0] # lifting condens. level TEMF
            self.thtemf  = np.apply_over_axes(np.mean,wrfin.variables["THUP_TEMF"][t0:t1,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # potential temperature updrafts [K]
            self.qttemf  = np.apply_over_axes(np.mean,wrfin.variables["QTUP_TEMF"][t0:t1,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # total water mixing rat. updrafts [kg kg-1]
            self.qltemf  = np.apply_over_axes(np.mean,wrfin.variables["QLUP_TEMF"][t0:t1,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # liq. water mixing rat. updrafts [K]
            self.c3dtemf = np.apply_over_axes(np.mean,wrfin.variables["CF3D_TEMF"][t0:t1,:,jj-n:jj+n1,ii-n:ii+n1],[2,3])[:,:,0,0] # cloud fraction updrafts [K]
            self.TEMF    = True
        except KeyError: # Other schemes don't have these variables, set to filval
            self.zi      = np.apply_over_axes(np.mean,wrfin.variables["PBLH"][t0:t1,jj-n:jj+n1,ii-n:ii+n1],[1,2])[:,0,0]
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
        self.slps          = self.ps / (1.-2.25577e-5 * self.hgt)**5.25588 # sea level pressure
        rhos               = self.ps / (Rd * self.T2) # surface density [kg m-3]
        wthvs              = (self.hfx/(rhos*cp)) + 0.61*self.T2*(self.lh/(rhos*Lv)) # surface buoyancy flux [W m-2]
        wthvs[wthvs<0]     = 0. # remove negative flux for w* calc
        self.zi[self.zi<0] = 0. # remove negative updraft heights (..)
        self.wstar         = (g * self.zi * wthvs/tref)**(1./3.) # convective velocity scale w*
        self.wstar[self.wstar < olga.sinkGlider] = 0. # convective velocity scale w* - sink glider

        #self.pV          = VgemCrossCountry(self.zi,self.wstar) # potential cross-country velocity
        #self.cPFD        = np.zeros_like(self.pV) # cumulative potential flight distance
        #for t in range(1,nt):
        #    self.cPFD[t] = self.cPFD[t-1] + av(self.pV[t-1],self.pV[t]) * self.dt
        #    if(self.hour[t] < self.hour[t-1]): # new day: reset PFD
        #        self.cPFD[t] = 0.

        # Get potential incoming shortwave radiation:
        self.swd_theory = np.zeros(nt)
        self.swd_frac = np.zeros(nt)
        for t in range(nt):
            self.swd_theory[t] = swin(self.doy[t],self.hour[t],self.lat,self.lon)
            self.swd_frac[t] = min([1.,max([1e-3,self.swd[t]]) / max([1.e-3,self.swd_theory[t]])])

        # Calculate temperature and dewpoint from potential temperature and total water mixing ratio
        self.qt  = self.qv + self.ql  # total water mixing ratio [kg kg-1]
        self.T   = np.zeros_like(self.th) # absolute temperature [K]
        self.Td  = np.zeros_like(self.th) # dew point temperature [K]
        self.Tu  = np.zeros_like(self.th) # absolute temperature updrafts [K]
        self.Tdu = np.zeros_like(self.th) # dew point temperature updrafts [K]

        e         = ((self.ps) * self.q2) / ((Rd/Rv) + self.q2) # vapor pressure
        self.Td2  = ((1./273.) - (Rv/Lv) * np.log(e/611.))**-1.

        for t in range(nt):
            self.th[t,:]  = self.th[t,:] + 300. # potential temp = base state (300) + perturbation (T)
            self.T[t,:]   = self.th[t,:] * (self.p[t,:] / 1.e5)**(287.05/1004.) # 1.e5 = reference pressure
            e             = ((self.p[t,:]) * self.qt[t,:]) / ((Rd/Rv) + self.qt[t,:]) # vapor pressure
            self.Td[t,:]  = ((1./273.) - (Rv/Lv) * np.log(e/611.))**-1.

            if(self.TEMF):
                self.Tu[t,:]  = self.thtemf[t,:] * (self.p[t,:] / 1.e5)**(287.05/1004.) # 1.e5 = reference pressure
                e             = ((self.p[t,:]) * self.qttemf[t,:]) / ((Rd/Rv) + self.qttemf[t,:]) # vapor pressure
                # ?!?!
                e[np.where(e<=0)] = 1e-9
                self.Tdu[t,:] = ((1./273.) - (Rv/Lv) * np.log(e/611.))**-1.
            else:
                self.Tu[t,:]  = -9999
                self.Tdu[t,:] = -9999


