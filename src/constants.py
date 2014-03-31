cp   = 1004.       # specific heat (constant P) [J kg-1 K-1]
g    = 9.8         # grav acceleration []
tref = 295.        # reference temperature [K]
rho  = 1.2         # density [kg m-3]
Rd   = 287.06      # gas constant dry air
Rv   = 461.5       # gas constant moits air
Lv   = 2.45e6      # Latent heat of vaporization
kappa = 0.4

#a    = 1.05e-4     # a-coef Ventus2 18, 40kg/m3
#b    = 1.64e-2     # b-coef Ventus2 18, 40kg/m3
#c    = 1.15        # c-coef Ventus2 18, 40kg/m3 

a    = 1.64e-4     # a-coef LS8-18, 40kg/m3
b    = 2.87e-2     # b-coef LS8-18, 40kg/m3
c    = 1.83        # c-coef LS8-18, 40kg/m3 

supd = 1.0         # Glider sink in updraft (..)
peff = 0.8        # efficiency pilot (..)

m2k = 1.95         # convert m/s to kts 

eps = 1.e-12       # small number 
