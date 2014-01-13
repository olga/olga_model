import numpy 
from pylab import *

close('all')

# Diagram settings
# ==============================
c1           = '#5aba39'   # Green
c2           = '#000000'   # Black
c4           = '#808080'   # Gray
c5           = '#93256d'   # Parcel 
a1           = 0.5         # Opacities
s1           = 9           # Font size 1
s2           = 10          # Font size 2
antia        = True        # Antialiasing?

# To-do: fix settings below, not working :(
p0           = 105000.     # highest pressure in diagram
p0_ref       = 100000.     # reference pressure 
ptop         = 10000.      # lowest pressue in diagram
top          = 100.        # ?
dp           = 100.        # ?
# ==============================

# Constants
# ------------------------------
Rd  = 287.
T0  = 273.15
e0  = 610.78
Lv  = 2.45e6 
Rv  = 461.
cp  = 1004.
g   = 9.81
eps = 0.622

# functions to skew and inverse skew diagram
# ------------------------------
def skewty(p):
  return (132.182 - 44.061 * log10(p))
def iskewty(y):
  return 10**((-y + 132.182)/44.061)
def skewtx(T,y):
  return (0.54 * T + 0.90692 * y)

# esat Tetens
# ------------------------------
def esat(Tin): 
  b    = 17.2694
  Tc   = Tin - T0 
  return e0 * numpy.exp((17.2693882 * Tc) / (Tc + 237.3))    

# esat Wobus
# ------------------------------
def esat_wob(Tin): 
  Tc   = Tin - T0
  e02  = e0 / 100.

  pol   = 0.99999683   + Tc*(-0.90826951E-02 + \
    Tc*(0.78736169E-04 + Tc*(-0.61117958E-06 +   \
    Tc*(0.43884187E-08 + Tc*(-0.29883885E-10 +   \
    Tc*(0.21874425E-12 + Tc*(-0.17892321E-14 +   \
    Tc*(0.11112018E-16 + Tc*(-0.30994571E-19)))))))))
  return (e02/pol**8.)

# qsat
# ------------------------------
def w(Tin,pin):
  x    = esat(Tin)
  return (eps * x) / (pin - x)

# Moist adiabatic lapse rate 
# ------------------------------
def dwobfskewt(Tin):
  xx = Tin - 20.
  if(xx <= 0.):
    pol = 1. + xx * (-8.8416605e-03 + \
               xx * (1.4714143e-04  + \
               xx * (-9.6719890e-07 + \
               xx * (-3.2607217e-08 + \
               xx * (-3.8598073e-10 )))))
    dwobfskewt = 15.130 / pol**4.
  else:
    pol = 1. + xx * (3.6182989e-03  + \
               xx * (-1.3603273e-05 + \
               xx * (4.9618922e-07  + \
               xx * (-6.1059365e-09 + \
               xx * (3.9401551e-11  + \
               xx * (-1.2588129e-13 + \
               xx * (1.6688280e-16)))))))
    dwobfskewt = 29.930 / pol**4 + 0.96 * xx - 14.8
  return dwobfskewt

def dsatlftskewt(tcin,pin):
  rdcp = Rd / cp

  if(pin == p0_ref):
    ma = tcin 
  else:
    pwrp  = (pin / p0_ref)**rdcp      
    tone  = ((tcin + T0) * pwrp) - T0
    eone  = dwobfskewt(tone) - dwobfskewt(tcin)

    rate  = 1.
    dlt   = 1.
    n     = 0
    while(abs(dlt) > 0.1):
      if(n != 0):
        rate = (ttwo - tone) / (etwo - eone)
        tone = ttwo
        eone = etwo
      ttwo  = tone - eone * rate
      pt    = (ttwo + T0) / pwrp - T0
      etwo  = pt + dwobfskewt(ttwo) - dwobfskewt(pt) - tcin
      dlt   = etwo * rate
      n     = n + 1
    ma = ttwo - dlt
  return ma

# Find closest value in list
# ------------------------------
def closest(target,collection):
  return min((abs(target - i), i) for i in collection)[1]

# Main routine to plot diagram
# ------------------------------
def skewtp(Tin,Tdin,pin,Tsin=[0],Tdsin=[0],psin=[0]):

  fig = figure(figsize=(8,9))
  #                   L    B    R    T    ws  hs
  fig.subplots_adjust(0.08,0.05,.99,0.95,0.2,0.08)
  subplot(111)

  # Create isotherms
  # ------------------------------
  isot = ([-100.,-90. ,-80. ,-70. ,-60. ,-50. ,-40. ,-30. , -20.,-10. ,  0. , 10. , 20. , 30.])
  lendt= ([ 132., 181., 247., 337., 459., 625., 855.,1050.,1050.,1050.,1050.,1050.,1050.,1050.])
  rendt= ([ top, top, top, top, top, top, top, 135., 185., 251., 342., 430., 556., 760.])

  for i in range(size(isot)):
    yy = ([skewty(lendt[i]),skewty(rendt[i])])
    xx = ([skewtx(isot[i],yy[0]),skewtx(isot[i],yy[1])])
    plot(xx,yy,color=c1,alpha=a1,aa=antia)
    if(i < 7):
      text(xx[1]+0.4,yy[1]+0.5,int(isot[i]),ha='center',color=c1,size=s1)
    else:
      text(xx[0]-0.5,yy[0]-1.5,int(isot[i]),ha='center',color=c2,size=s2)

  # Create isobars
  # ------------------------------
  pres = ([1050.,1000.,850.,700.,500.,400.,300.,250.,200.,150.,100.])
  xpl  = ([-19.0,-19.0,-19.0,-19.0,-19.0,-19.0,-19.0,-19.0,-19.0,-19.0,-19.0])
  xpr  = ([20.1,20.1,20.1,20.1,20.83,18.6,18.6,18.6,18.6,18.6,18.6])

  for i in range(1,size(pres)):
    xx = ([xpl[i],xpr[i]])
    yy = ([skewty(pres[i]),skewty(pres[i])])
    plot(xx,yy,color=c1,alpha=a1,aa=antia)
    text(xx[0]-0.5,yy[0],int(pres[i]),ha='right',va='center',color=c2,size=s2)

  # Create dry adiabats
  # ------------------------------
  thetas= ([-30. ,-20. ,-10. ,0.   ,10.  ,20.  ,30.  ,40. ,50. ,60. ,70. ,80.]) 
  lendth= ([880. ,670. ,512. ,388. ,292. ,220. ,163. ,119.,100.,100.,100.,100.])
  rendth= ([1000.,1000.,1000.,1000.,1000.,1000.,1000.,980.,825.,695.,590.,500.])
  
  dp1 = 1.
  for i in range(size(thetas)):
    np      = int((rendth[i]-lendth[i])/dp1)	    
    theta   = numpy.zeros(np)
    ps      = numpy.zeros(np)

    p1      = rendth[i]
    for j in range(np):
      ps[j]    = skewty(p1)      
      theta[j] = skewtx(((thetas[i]+T0) * (p1 / (p0_ref/100.))**(Rd/cp) - T0),ps[j])
      p1       = p1 - dp1
    plot(theta,ps,color=c4)
    text(theta[-1],ps[-1]-1.,int(thetas[i]),color=c4,backgroundcolor='w',size=s1,ha='left')

  # Create moist adiabats
  # ------------------------------
  pseudo   = ([32.,28.,24.,20.,16.,12.,8.])
  p1        = arange(p0_ref,23000,-100)
  np       = size(p1)
  ps       = numpy.zeros(np)
  ma       = numpy.zeros(np)

  for i in range(size(pseudo)):
    for j in range(np):
      moista = dsatlftskewt(pseudo[i],p1[j])
      ps[j]  = skewty(p1[j]/100.)
      ma[j]  = skewtx(moista,ps[j])

    plot(ma,ps,'--',color=c4)
    text(ma[-1]-0.7,ps[-1],int(pseudo[i]),color=c4,backgroundcolor='w',size=s1,ha='left')

  # Create isohumes / mixing ratio
  # ------------------------------
  mixrat   = ([20.,12.,8.,5.,3.,2.,1.])
  p1        = arange(p0_ref,40000,-100)
  np       = size(p1)
  ps       = numpy.zeros(np)
  mr       = numpy.zeros(np)

  for i in range(size(mixrat)):
    for j in range(np):
      mix  = (((1. / T0) - (Rv / Lv) * log(((mixrat[i]/1000.) * (p1[j]/1000.)) / (0.611 * ((mixrat[i]/1000.)+0.622))))**-1.) - T0
      ps[j] = skewty(p1[j]/100.)
      mr[j] = skewtx(mix,ps[j])
    plot(mr,ps,'k:')
    text(mr[-1]-0.0,ps[-1],int(mixrat[i]),backgroundcolor='w',size='9',ha='center')

  # Plot sounding measurement
  # ------------------------------
  T_skew  = numpy.zeros_like(Tin)
  Td_skew = numpy.zeros_like(Tin)
  p_skew  = numpy.zeros_like(pin)

  if(Tsin.size>1):
    Ts_skew = numpy.zeros_like(Tsin)
    Tds_skew = numpy.zeros_like(Tsin)
    ps_skew = numpy.zeros_like(Tsin)

    for i in range(Ts_skew.size):
      ps_skew[i]  = skewty(psin[i])
      Ts_skew[i]  = skewtx(Tsin[i],ps_skew[i])
      Tds_skew[i] = skewtx(Tdsin[i],ps_skew[i])

  for i in range(p_skew.size):
    p_skew[i]  = skewty(pin[i]/100.)
    T_skew[i]  = skewtx(Tin[i]-T0,p_skew[i])
    Td_skew[i] = skewtx(Tdin[i]-T0,p_skew[i])



  plot(T_skew,p_skew,'r-',linewidth=2)
  plot(Td_skew,p_skew,'b-',linewidth=2)

  if(Tsin.size>1):
    plot(Ts_skew,ps_skew,'r--',linewidth=1.5,dashes=[5,2])
    plot(Tds_skew,ps_skew,'b--',linewidth=1.5,dashes=[5,2])


  # _________________________________________________________
  # Plot atmospheric data
  #plot(T_skew,p_skew,color='r',linewidth=2,aa=antia)
  #plot(Td_skew,p_skew,color='b',linewidth=2,aa=antia)

  ## Wind barbs
  #x_barb = 26.
  #plot(([x_barb,x_barb]),([skewty(p[0]/100.),skewty(p[-1]/100.)]),color=c4)
  #for i in range(size(p)):
  #  if((i < 30) & (i%4 == 0.)):	    
  #    barbs(x_barb,skewty(p[i]/100.),u[i],v[i],color=c4,length=6)
  #    text(x_barb-2.5,skewty(p[i]/100.),int(z[i]),color=c4,size=s1,va='center',ha='right')
  #  elif(i>30):    
  #    barbs(x_barb,skewty(p[i]/100.),u[i],v[i],color=c4,length=6)
  #    text(x_barb-2.5,skewty(p[i]/100.),int(z[i]),color=c4,size=s1,va='center',ha='right')

  ##y_line = skewty(p[-1]/100.)+1.
  ##plot(([x_barb-5.0,x_barb]),([y_line,y_line]),'-',color=c4)
  ##text(x_barb-2.5,y_line+0.5,'(m)  |  (kts)',color=c4,size=s1,ha='center')

  ## Parcels!
  ## Surface based parcel
  #par = parcel(T[0],q[0],p[0],0)
  #plot(par[5,0:par[6,0]],par[3,0:par[6,0]],'-',color=c5,linewidth=1)
  #plot(par[4,:],par[3,:],'-',color=c5,linewidth=1)
  #fill_betweenx(par[3,:],par[8,:],par[4,:],where=par[4,:]>=par[8,:],color=c1,alpha=0.3)

  #info = ([
  #        r'$CAPE_{SB} = $' + '%.0f' % par[6,2],
  #        r'$LI_{SB} = $' + '%.2f' % par[6,1],
  #        ])
  #xt = 18.6
  #yt = 42.75

  #for k in range(size(info)):
  #  text(xt,yt,info[k],ha='right',size=s1,color=c2)
  #  yt = yt - 1

  # _________________________________________________________
  # Remove frame, labels, draw manually
  xticks([])
  yticks([])
  box('off')
  xmin = skewtx(-109.1,skewty(100.))
  xmax = skewtx(51.75, skewty(1050.)) -6.
  ymin = skewty(1050.)
  ymax = skewty( 100.)

  xlim(xmin-1,xmax+8)
  ylim(ymin-1,ymax+2)
  
  plot(([xmin,xmin]),([ymin,ymax]),'k-',aa=antia)
  plot(([xmin,xmax]),([ymin,ymin]),'k-',aa=antia)

  figtext(0.5,0.018,'Temperature (C)',size=s2,ha='center')
  figtext(0.025,0.52,'Pressure (hPa)',rotation=90,size=s2)
  label = 'Skew-T log-P, LOCATION, XXXXUTC' 
  figtext(0.5,0.97,label,size=11,ha='center')
  label = '18 x 18 km GFS-initiated WRF-ARW'
  figtext(0.5,0.95,label,size=8,ha='center')

  #savefig('skewtp.png')



#skewtp()


