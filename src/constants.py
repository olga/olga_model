#
# Copyright (c) 2013-2016 Bart van Stratum
# Copyright (c) 2015-2016 Roel Baardman
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

cp   = 1004.       # specific heat (constant P) [J kg-1 K-1]
g    = 9.8         # grav acceleration []
tref = 295.        # reference temperature [K]
rho  = 1.2         # density [kg m-3]
Rd   = 287.06      # gas constant dry air
Rv   = 461.5       # gas constant moits air
Lv   = 2.45e6      # Latent heat of vaporization
kappa = 0.4        # Von Karman constant
m2k = 1.95         # convert m/s to kts 
eps = 1.e-12       # small number 
filval = -9999     # Fill value for missing data
