import numpy as np
from pylab import *

close('all')

def Vgem(wg, a, b, c):
    vstf    = -(b-(b-(4.*-a*(-c-wg))**0.5))/(2.*-a) # speed to fly given updraft (MacCready) velocity [km h-1]
    wstf    = -a*vstf**2.+b*vstf-c # sink glider at Vstf [m s-1]
    alpha   = -wstf / (wg - wstf) # fraction time spent circling [-]
    Vgem    = (1.-alpha) * vstf # (1-alpha)*V = cross-country speed [km h-1] 
    return Vgem

# Polar coefficients: ventus2-18 at 45 kg/m3, LS8-15 at 45 kg/m3, std Cirrus at 30 kg/m3
a  = ([1.05e-4, 2.14e-4, 1.71e-4]) 
b  = ([1.79e-2, 4.65e-2, 2.43e-2])
c  = ([1.3099,  3.2768,  1.4600 ])
t0 = 8.   # Start convection
t1 = 18.  # End convection
dt = 0.1  # Time step calculation

time  = np.arange(t0, t1+0.01, dt)
wstar = 5*np.sin(np.pi*(time-t0)/(t1-t0)) # Updraft velocity
wglider = wstar - 0.7 # Assume 1m/s sink glider in updraft
wglider[wglider<0] = 0.

v1 = Vgem(wglider, a[0], b[0], c[0]) # Ventus2 
v2 = Vgem(wglider, a[1], b[1], c[1]) # Std Cirrus
v3 = Vgem(wglider, a[2], b[2], c[2]) # Std Cirrus

pfd1 = np.zeros_like(v1)
pfd2 = np.zeros_like(v2)
pfd3 = np.zeros_like(v3)

for t in range(1,pfd1.size):
    pfd1[t] = pfd1[t-1] + v1[t] * dt 
    pfd2[t] = pfd2[t-1] + v2[t] * dt 
    pfd3[t] = pfd3[t-1] + v3[t] * dt 

# Simple check:
print('PFD Ventus = %f km, LS8 = %f km, Cirrus = %f km'%(np.sum(v1)*dt, np.sum(v2)*dt, np.sum(v3)*dt))

figure()
subplot(131)
plot(time,wstar,label='updraft velocity')
plot(time,wglider,label='vertical velocity glider in updraft')
xlabel('time UTC')
ylabel('Vertical velocity w [m/s]')
legend(frameon=False)

subplot(132)
plot(time,v1,label='Ventus2-18 at 45 kg/m3')
plot(time,v2,label='LS8-15 at 45 kg/m3')
plot(time,v3,label='Std-Cirrus at 30 kg/m3')
xlabel('time UTC')
ylabel('cross-country velocity [km/h]')
legend(frameon=False)

subplot(133)
plot(time,pfd1,label='Ventus2-18 at 45 kg/m3')
plot(time,pfd2,label='LS8-15 at 45 kg/m3')
plot(time,pfd3,label='Std-Cirrus at 30 kg/m3')
xlabel('time UTC')
ylabel('Cumulative potential flight distance [km]')
legend(frameon=False)
