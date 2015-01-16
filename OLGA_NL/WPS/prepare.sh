#!/bin/bash 
# Script to prepare directory for WPS routines

# -------------------------
# Change WPSPath to relative or absolute path to WPS installation directory
#WPSPath="../../../WPSv361"  # Mint
WPSPath="../../../WPS_v3.6.1" # Thunder
# -------------------------

# Link the required files to current directory
ln -s $WPSPath/*.exe .
ln -s $WPSPath/link_grib.csh .
ln -s $WPSPath/geogrid .
ln -s $WPSPath/metgrid .
ln -s $WPSPath/ungrib/Variable_Tables/Vtable.GFS Vtable
