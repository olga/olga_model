import numpy as np
from netCDF4 import Dataset 
import sys
import os

if(len(sys.argv)<2):
  sys.exit('specify input wrfout file')
else:
  filein  = sys.argv[1]
  print('processing file %s'%filein)
  wrfin   = Dataset(filein,'r')
  fileout = filein+'_red.nc'
  wrfout  = Dataset(fileout,'w',format='NETCDF4')

  # Dimensions input file
  nt    = len(wrfin.variables["XTIME"][:])
  we    = len(wrfin.dimensions["west_east"])  
  wes   = len(wrfin.dimensions["west_east_stag"])  
  sn    = len(wrfin.dimensions["south_north"])  
  sns   = len(wrfin.dimensions["south_north_stag"])  
  bt    = len(wrfin.dimensions["bottom_top"])
  bts   = len(wrfin.dimensions["bottom_top_stag"])
  sls   = len(wrfin.dimensions["soil_layers_stag"])
  dstrl = len(wrfin.dimensions["DateStrLen"])

  # Create dimensions in output file 
  Time             = wrfout.createDimension('Time',nt)
  west_east        = wrfout.createDimension('west_east',we)
  west_east_stag   = wrfout.createDimension('west_east_stag',wes)
  south_north      = wrfout.createDimension('south_north',sn)
  south_north_stag = wrfout.createDimension('south_north_stag',sns)
  bottom_top       = wrfout.createDimension('bottom_top',bt)
  bottom_top_stag  = wrfout.createDimension('bottom_top_stag',bts)
  soil_layers_stag = wrfout.createDimension('soil_layers_stag',sls)
  DateStrLen       = wrfout.createDimension('DateStrLen',dstrl)

  # Define variables
  T00              = wrfout.createVariable('T00','f4',('Time',))
  P00              = wrfout.createVariable('P00','f4',('Time',))
  PH               = wrfout.createVariable('PH','f4',('Time','bottom_top_stag','south_north','west_east',))
  PHB              = wrfout.createVariable('PHB','f4',('Time','bottom_top_stag','south_north','west_east',))

  T2               = wrfout.createVariable('T2','f4',('Time','south_north','west_east'))
  U10              = wrfout.createVariable('U10','f4',('Time','south_north','west_east'))
  V10              = wrfout.createVariable('V10','f4',('Time','south_north','west_east'))
  PSFC             = wrfout.createVariable('PSFC','f4',('Time','south_north','west_east'))
  RAINNC           = wrfout.createVariable('RAINNC','f4',('Time','south_north','west_east'))
  RAINC            = wrfout.createVariable('RAINC','f4',('Time','south_north','west_east'))
  RAINSH           = wrfout.createVariable('RAINSH','f4',('Time','south_north','west_east'))

  PBLH             = wrfout.createVariable('PBLH','f4',('Time','south_north','west_east'))
  LANDMASK         = wrfout.createVariable('LANDMASK','f4',('Time','south_north','west_east')) 
  HFX              = wrfout.createVariable('HFX','f4',('Time','south_north','west_east')) 
  LH               = wrfout.createVariable('LH','f4',('Time','south_north','west_east')) 

  CLDFRA           = wrfout.createVariable('CLDFRA','f4',('Time','bottom_top','south_north','west_east'))
  T                = wrfout.createVariable('T','f4',('Time','bottom_top','south_north','west_east'))
  #W                = wrfout.createVariable('W','f4',('Time','bottom_top_stag','south_north','west_east'))

  Times            = wrfout.createVariable('Times','S1',('Time','DateStrLen'))
 
  try:
    dummy = wrfin.variables["WUPD_TEMF"][0,0,0,0]
    temf  = True
  except KeyError:
    print 'no TEMF'
    temf  = False

  if(temf):
    WUPD_TEMF      = wrfout.createVariable('WUPD_TEMF','f4',('Time','bottom_top','south_north','west_east'))
    CF3D_TEMF      = wrfout.createVariable('CF3D_TEMF','f4',('Time','bottom_top','south_north','west_east'))
    HD_TEMF        = wrfout.createVariable('HD_TEMF','f4',('Time','south_north','west_east'))
    HCT_TEMF       = wrfout.createVariable('HCT_TEMF','f4',('Time','south_north','west_east'))
    CFM_TEMF       = wrfout.createVariable('CFM_TEMF','f4',('Time','south_north','west_east'))
    LCL_TEMF       = wrfout.createVariable('LCL_TEMF','f4',('Time','south_north','west_east'))

  # Copy variables
  T00[:]           = wrfin.variables["T00"][:] 
  P00[:]           = wrfin.variables["P00"][:]
  PH[:,:,:,:]      = wrfin.variables["PH"][:,:,:,:]
  PHB[:,:,:,:]     = wrfin.variables["PHB"][:,:,:,:]
                     
  T2[:,:,:]        = wrfin.variables["T2"][:,:,:]
  U10[:,:,:]       = wrfin.variables["U10"][:,:,:]
  V10[:,:,:]       = wrfin.variables["V10"][:,:,:]
  PSFC[:,:,:]      = wrfin.variables["PSFC"][:,:,:]
  RAINNC[:,:,:]    = wrfin.variables["RAINNC"][:,:,:]
  RAINC[:,:,:]     = wrfin.variables["RAINC"][:,:,:]
  RAINSH[:,:,:]    = wrfin.variables["RAINSH"][:,:,:]
                     
  PBLH[:,:,:]      = wrfin.variables["PBLH"][:,:,:]
  LANDMASK[:,:,:]  = wrfin.variables["LANDMASK"][:,:,:]
  HFX[:,:,:]       = wrfin.variables["HFX"][:,:,:]
  LH[:,:,:]        = wrfin.variables["LH"][:,:,:]
                     
  CLDFRA[:,:,:,:]  = wrfin.variables["CLDFRA"][:,:,:,:]
  T[:,:,:,:]       = wrfin.variables["T"][:,:,:,:]
  #W[:,:,:,:]       = wrfin.variables["W"][:,:,:,:]
                     
  Times[:,:]       = wrfin.variables["Times"][:,:]

  if(temf):
    WUPD_TEMF[:,:,:,:] = wrfin.variables["WUPD_TEMF"][:,:,:,:]
    CF3D_TEMF[:,:,:,:] = wrfin.variables["CF3D_TEMF"][:,:,:,:]
    HD_TEMF[:,:,:]     = wrfin.variables["HD_TEMF"][:,:,:]
    HCT_TEMF[:,:,:]    = wrfin.variables["HCT_TEMF"][:,:,:]
    CFM_TEMF[:,:,:]    = wrfin.variables["CFM_TEMF"][:,:,:]
    LCL_TEMF[:,:,:]    = wrfin.variables["LCL_TEMF"][:,:,:]

  wrfin.close()
  wrfout.close()

  print 'size old-new (MiB) =',os.path.getsize(filein)/(1024.*1024.),os.path.getsize(fileout)/(1024.*1024.)
