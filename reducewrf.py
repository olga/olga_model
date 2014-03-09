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
  fileout = filein+'_2d.nc'
  wrfout  = Dataset(fileout,'w',format='NETCDF4')

  # Dimensions input file
  nt                = len(wrfin.variables["XTIME"][:])
  we                = len(wrfin.dimensions["west_east"])  
  wes               = len(wrfin.dimensions["west_east_stag"])  
  sn                = len(wrfin.dimensions["south_north"])  
  sns               = len(wrfin.dimensions["south_north_stag"])  
  bt                = len(wrfin.dimensions["bottom_top"])
  bts               = len(wrfin.dimensions["bottom_top_stag"])
  sls               = len(wrfin.dimensions["soil_layers_stag"])
  dstrl             = len(wrfin.dimensions["DateStrLen"])

  # Create dimensions in output file 
  Time              = wrfout.createDimension('Time',nt)
  west_east         = wrfout.createDimension('west_east',we)
  west_east_stag    = wrfout.createDimension('west_east_stag',wes)
  south_north       = wrfout.createDimension('south_north',sn)
  south_north_stag  = wrfout.createDimension('south_north_stag',sns)
  bottom_top        = wrfout.createDimension('bottom_top',bt)
  bottom_top_stag   = wrfout.createDimension('bottom_top_stag',bts)
  bottom_top_cloud  = wrfout.createDimension('bottom_top_cloud',3)
  soil_layers_stag  = wrfout.createDimension('soil_layers_stag',sls)
  DateStrLen        = wrfout.createDimension('DateStrLen',dstrl)

  # Define variables
  T00               = wrfout.createVariable('T00','f4',('Time',))
  P00               = wrfout.createVariable('P00','f4',('Time',))
  T2                = wrfout.createVariable('T2','f4',('Time','south_north','west_east'))
  U10               = wrfout.createVariable('U10','f4',('Time','south_north','west_east'))
  V10               = wrfout.createVariable('V10','f4',('Time','south_north','west_east'))
  PSFC              = wrfout.createVariable('PSFC','f4',('Time','south_north','west_east'))
  RAINNC            = wrfout.createVariable('RAINNC','f4',('Time','south_north','west_east'))
  RAINC             = wrfout.createVariable('RAINC','f4',('Time','south_north','west_east'))
  LANDMASK          = wrfout.createVariable('LANDMASK','f4',('south_north','west_east')) 
  HFX               = wrfout.createVariable('HFX','f4',('Time','south_north','west_east')) 
  LH                = wrfout.createVariable('LH','f4',('Time','south_north','west_east')) 
  Times             = wrfout.createVariable('Times','S1',('Time','DateStrLen'))
  # Derived
  CLDFRA            = wrfout.createVariable('CLDFRA','f4',('Time','bottom_top_cloud','south_north','west_east'))

  # See if TEMF available (small domain) 
  try:
    dummy = wrfin.variables["WUPD_TEMF"][0,0,0,0]
    temf  = True
  except KeyError:
    print 'no TEMF'
    temf  = False

  if(temf):
    HD_TEMF         = wrfout.createVariable('HD_TEMF','f4',('Time','south_north','west_east')) 
    HCT_TEMF        = wrfout.createVariable('HCT_TEMF','f4',('Time','south_north','west_east'))
    LCL_TEMF        = wrfout.createVariable('LCL_TEMF','f4',('Time','south_north','west_east'))
  else:
    PBLH            = wrfout.createVariable('PBLH','f4',('Time','south_north','west_east'))

  # Copy variables
  T00[:]            = wrfin.variables["T00"][:] 
  P00[:]            = wrfin.variables["P00"][:]
  T2[:,:,:]         = wrfin.variables["T2"][:,:,:]
  U10[:,:,:]        = wrfin.variables["U10"][:,:,:]
  V10[:,:,:]        = wrfin.variables["V10"][:,:,:]
  PSFC[:,:,:]       = wrfin.variables["PSFC"][:,:,:]
  RAINNC[:,:,:]     = wrfin.variables["RAINNC"][:,:,:]
  RAINC[:,:,:]      = wrfin.variables["RAINC"][:,:,:]
  LANDMASK[:,:]     = wrfin.variables["LANDMASK"][0,:,:]
  HFX[:,:,:]        = wrfin.variables["HFX"][:,:,:]
  LH[:,:,:]         = wrfin.variables["LH"][:,:,:]
  Times[:,:]        = wrfin.variables["Times"][:,:]

  if(temf):
    HD_TEMF[:,:,:]  = wrfin.variables["HD_TEMF"][:,:,:]
    HCT_TEMF[:,:,:] = wrfin.variables["HCT_TEMF"][:,:,:]
    LCL_TEMF[:,:,:] = wrfin.variables["LCL_TEMF"][:,:,:]
  else:
    PBLH[:,:,:]     = wrfin.variables["PBLH"][:,:,:]

  CLDFRA[:,0,:,:]     = np.sum(wrfin.variables["CLDFRA"][:,:35,:,:],axis=1)   / 35. 
  CLDFRA[:,1,:,:]     = np.sum(wrfin.variables["CLDFRA"][:,35:52,:,:],axis=1) / 17. 
  CLDFRA[:,2,:,:]     = np.sum(wrfin.variables["CLDFRA"][:,52:,:,:],axis=1)   / 11.

  wrfin.close()
  wrfout.close()

  print 'size old-new (MiB) =',os.path.getsize(filein)/(1024.*1024.),os.path.getsize(fileout)/(1024.*1024.)
