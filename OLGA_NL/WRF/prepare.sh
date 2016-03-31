#!/bin/bash 
# Script to prepare directory for WRF routines

# -------------------------
# Change WPSPath to relative or absolute path to WPS installation directory
WRFPath="/home/zmaw/m300241/WRFnl/WRFv361" # MPIPC
#WRFPath="../../../WRFv361" # Mint
#WRFPath="../../../WRF_v3.6.1" # Thunder
# -------------------------

# Cleanup
rm -rf wrf.exe
rm -rf real.exe
rm -rf ndown.exe
rm -rf nup.exe
rm -rf tc.exe
rm -rf LANDUSE.TBL
rm -rf VEGPARM.TBL
rm -rf SOILPARM.TBL
rm -rf GENPARM.TBL
rm -rf RRTM*
rm -rf CLM*
rm -rf *.formatted

# Link the required files to current directory
ln -s $WRFPath/main/wrf.exe .
ln -s $WRFPath/main/real.exe .
ln -s $WRFPath/main/ndown.exe .
ln -s $WRFPath/main/nup.exe .
ln -s $WRFPath/main/tc.exe .
ln -s $WRFPath/run/LANDUSE.TBL .
ln -s $WRFPath/run/VEGPARM.TBL .
ln -s $WRFPath/run/SOILPARM.TBL .
ln -s $WRFPath/run/GENPARM.TBL .
ln -s $WRFPath/run/RRTM* .
ln -s $WRFPath/run/CLM* .
ln -s $WRFPath/run/*.formatted .
