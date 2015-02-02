#!/bin/bash 
# Script to prepare directory for WPS routines

# -------------------------
# Change WPSPath to relative or absolute path to WPS installation directory
WPSPath="/home/zmaw/m300241/WRFnl/WPSv361"  # MPIPC
#WPSPath="../../../WPSv361"  # Mint
#WPSPath="../../../WPS_v3.6.1" # Thunder
# -------------------------

# Cleanup
rm -rf geogrid
rm -rf geogrid.exe
rm -rf metgrid
rm -rf metgrid.exe
rm -rf ungrib
rm -rf ungrib.exe
rm -rf link_grib.csh

# Link the required files to current directory
ln -s $WPSPath/geogrid .
ln -s $WPSPath/geogrid.exe .
ln -s $WPSPath/metgrid .
ln -s $WPSPath/metgrid.exe .
ln -s $WPSPath/ungrib .
ln -s $WPSPath/ungrib.exe .
ln -s $WPSPath/link_grib.csh


