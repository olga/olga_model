#!/bin/bash 
# Script to prepare directory for WRF routines

# -------------------------
# Change WPSPath to relative or absolute path to WPS installation directory
WRFPath="../../../WRFv361" # Mint
#WRFPath="../../../WRF_v3.6.1" # Thunder
# -------------------------

# Link the required files to current directory
ln -s $WRFPath/main/*.exe .
ln -s $WRFPath/run/LANDUSE.TBL .
ln -s $WRFPath/run/VEGPARM.TBL .
ln -s $WRFPath/run/SOILPARM.TBL .
ln -s $WRFPath/run/GENPARM.TBL .
ln -s $WRFPath/run/RRTM* .
ln -s $WRFPath/run/CLM* .
