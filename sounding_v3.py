"""
Plot skewT-log(p) diagrams 
Bart van Stratum, 03-2014
"""

import numpy as np
import pylab as pl
import sys
import matplotlib.image as image

from matplotlib import rc
rc('font', size=9)
rc('legend', fontsize=8)

"""
Sounding settings / colors
"""
c1 = '#5aba39'        # Green
c3 = '#ff7f00'        # Orange
c5 = '#377eb8'        # Blue

c2 = '#000000'        # Black
c4 = '#808080'        # Gray

cT = '#e41a1c'        # red - T 
cTd = '#386cb0'       # blue - Td  

a1 = 0.2              # Opacity
a2 = 0.5              # Opacity

"""
Constants
"""
Rd  = 287.     
Rv  = 461.     
eps = Rd/Rv
Lv  = 2.45e6   
cp  = 1004.    
T0  = 273.15    
p0  = 1.e5    
e0  = 0.611 
g   = 9.81

"""
Given pressure (pa), returns y in figure coordinates (plus function for its inverse)
"""
def skewty(pin):
    return (132.182 - 44.061 * np.log10(pin))
def iskewty(yin):
    return 10**((-yin + 132.182)/44.061)

"""
Given Temperature (K) and y (fig coords, output skewty), returns x in fig coords
"""
def skewtx(Tin,yin):
    return (0.54 * Tin + 0.90692 * yin)
def iskewtxy(sk,yin):
    return (sk - 0.90692 * yin) / 0.54
def iskewtxT(sk,Tin):
    return (sk - 0.54 * Tin) / 0.90692

"""
Calculate saturation vapor pressure (Wobus)
"""
def es(Tin): 
    Tc = Tin - T0
    e02 = e0 * 10.
    pol = 0.99999683   + Tc*(-0.90826951E-02 + \
     Tc*(0.78736169E-04 + Tc*(-0.61117958E-06 + \
     Tc*(0.43884187E-08 + Tc*(-0.29883885E-10 + \
     Tc*(0.21874425E-12 + Tc*(-0.17892321E-14 + \
     Tc*(0.11112018E-16 + Tc*(-0.30994571E-19)))))))))
    return (e02/pol**8.)

"""
Calculate saturation mixing ratio from T and p
"""
def rs(Tin,pin):
    x = es(Tin) * 100.
    return (eps * x) / (pin - x)

"""
Calculate vapor pressure from p and mixing ratio 
"""
def r2es(rin,pin):
    return (pin * rin) / (eps + rin)

""" 
Calculate dew point temperature
"""
def Td(rin,pin):
    e = r2es(rin,pin)
    return ((1./T0)-(Rv/Lv)*np.log(e/e0))**-1.  

"""
Moist adiabatic lapse rate
"""
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

"""
Moist adiabatic lapse rate
"""
def dsatlftskewt(tcin,pin):
    pwrp = (pin / p0)**(Rd/cp)      
    tone = ((tcin + T0) * pwrp) - T0
    eone = dwobfskewt(tone) - dwobfskewt(tcin)

    rate = 1.
    dlt = 1.
    n = 0
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

"""
Find closest value in list
"""
def closest(target,collection):
    return min((abs(target - i), i) for i in collection)[1]

"""
Exner function
"""
def exner(pin,pref=p0):
    return (pin/pref)**(Rd/cp)

"""
Input class for sounding
"""
class skewt_input:
    def __init__(self):
        self.stype = 0          # 0=high top, 1=low top
        self.T  = np.array([])  # Sounding temperature
        self.Td = np.array([])  # Sounding dewpoint temperature
        self.u  = np.array([])  # Sounding u-wind
        self.v  = np.array([])  # Sounding v-wind
        self.p  = np.array([])  # Sounding pressure levels
        self.z  = np.array([])  # Sounding height levels
        self.name = ""          # Location
        self.time = ""          # Time
        self.parcel = False     # Lauch parcel?
        self.ps = -1            # Surface pressure
        self.Ts = -9999         # Surface temperature of parcel
        self.rs = -9999         # Surface moisture mixing ratio of parcel

"""
Main routine which makes the diagram
input: si = skewt_input object
"""
def skewtlogp(si):
    """
    Settings for high and low sounding top
    """
    if(si.stype==0):
        pbottom = 105000.       # highest pressure in diagram (bottom)
        ptop = 10000.           # lowest pressue in diagram (top)
        ptop_thd = 40000        # top of dry adiabats
        ptop_mxr = 60000        # top of mixing ratio lines
        pbottom_thw = pbottom   # bottom of moist adiabats
        dp = 100.               # pressure interval used in some calculations
        Tleft = -35.            # lowest temperature (@pb) in diagram (left)
        Tright = 35.            # highest temperature (@pb) in diagram (right)
        dp_label = 7000

        isotherms     = np.arange(-100,30.001,10)
        isobars       = np.array([1050.,1000.,850.,700.,500.,400.,300.,250.,200.,150.,100.,50])*100.
        dryadiabats   = np.arange(-30,50.001,10) 
        moistadiabats = np.array([28.,24.,20.,16.,12.,8.])
        mixrat        = np.array([20.,12.,8.,5.,3.,2.,1.])
    elif(si.stype==1):
        pbottom = 105000.       # highest pressure in diagram (bottom)
        ptop = 50000.           # lowest pressue in diagram (top)
        ptop_thd = 70000        # top of dry adiabats
        ptop_mxr = 70000        # top of mixing ratio lines
        pbottom_thw = pbottom   # bottom of moist adiabats
        dp = 100.               # pressure interval used in some calculations
        Tleft = -11.            # lowest temperature (@pb) in diagram (left)
        Tright = 35.            # highest temperature (@pb) in diagram (right)
        dp_label = 2500         # spacing (in pressure coords) of some label placement

        isotherms     = np.arange(-100,60.001,10)
        isobars       = np.array([1050.,1000.,900.,850.,700.,600.,500.])*100.
        dryadiabats   = np.arange(-10,30.001,5) 
        moistadiabats = np.array([28.,24.,20.,16.,12.,8.])
        mixrat        = np.array([20.,12.,8.,5.,3.,2.,1.])
    else:
        sys.exit('stype=%i not supported!'%stype)

    """
    Setup figure
    """
    fig = pl.figure(figsize=(6.5,6))
    #                   L    B    R    T    ws  hs
    fig.subplots_adjust(0.08,0.10,1.0,0.93,0.2,0.08)
    pl.subplot(111)

    # Calculate bounds in figure coordinates
    y00 = skewty(pbottom)
    y11 = skewty(ptop)
    x00 = skewtx(Tleft,y00)
    x11 = skewtx(Tright,y00)

    # Spacing for some labels
    hs = np.abs(x11-x00)/200.
    vs = np.abs(y11-y00)/200.
 
    """
    1. Create isotherms
    """
    for T in isotherms:
        y = ([skewty(pbottom),skewty(ptop)])
        x = ([skewtx(T,y[0]),skewtx(T,y[1])])

        # Check if partially out of bounds
        lloc = 0
        if(x[0] < x00):
            x[0] = x00
            y[0] = iskewtxT(x00,T)
        if(x[1] > x11):
            x[1] = x11
            y[1] = iskewtxT(x11,T)

        if(x[1] >= x00 and y[1] >= y00):
            pl.plot(x,y,color=c2,alpha=a1)
            if(x[0] > x00 and x[0] < x11):
              pl.text(x[0],y[0]-2*hs,int(T),ha='center',va='top',color=c2)
            else:
              pl.text(x[1]-5*hs,y[1]-5*vs,int(T),ha='center',va='center',color=c2,alpha=a2)

    """
    2. Create isobars
    """
    for p in isobars:
        # Check if out of bounds
        if((p >= ptop) and (p <= pbottom)):
            y = skewty(p)
            pl.plot([x00,x11],[y,y],color=c2,alpha=a1)
            pl.text(x00-hs,y,int(p/100.),ha='right',va='center',color=c2)
         
    """
    3. Create dry adiabats
    """
    p = np.arange(pbottom,(np.max(([ptop,ptop_thd]))-dp),-dp)
    x = np.zeros(p.size)   
    y = np.zeros(p.size)   
   
    for k in range(p.size):
        y[k] = skewty(p[k])
    for th in dryadiabats: 
        x[:] = 0.
        for k in range(p.size):
            xtmp = skewtx(((th+T0) * exner(p[k]))-T0,y[k])
            if(xtmp >= x00 and xtmp <= x11):
                x[k] = xtmp
            else:
                x[k] = -9999
    
        doPlot = np.where(x != -9999)
        if(len(doPlot[0]) > 0):
            lloc = max(0,np.size(x[doPlot])-int(5000./dp))
            pl.plot(x[doPlot],y[doPlot],color=c1)
            pl.text(x[doPlot][lloc],y[doPlot][lloc],int(th),color=c1,ha='center',va='center',backgroundcolor='w')

    """ 
    4. Create moist adiabats
    """
    p = np.arange(np.min(([pbottom,pbottom_thw])),ptop-dp,-dp)
    x = np.zeros(p.size)
    y = np.zeros(p.size)

    for thw0 in moistadiabats:
        for k in range(p.size):
            thw = dsatlftskewt(thw0,p[k])
            y[k] = skewty(p[k])
            xtmp = skewtx(thw,y[k])
            if(xtmp >= x00 and xtmp <= x11):
                x[k] = xtmp
            else:
                x[k] = -9999
    
        doPlot = np.where(x != -9999)
        if(len(doPlot[0]) > 0):
            lloc = max(0,np.size(x[doPlot])-int(8000./dp))
            pl.plot(x[doPlot],y[doPlot],'--',color=c3)
            pl.text(x[doPlot][lloc],y[doPlot][lloc],int(thw0),color=c3,ha='center',va='center',backgroundcolor='w')
 
    """ 
    5. Create isohumes / mixing ratio lines
    """
    p = np.arange(pbottom,(np.max(([ptop,ptop_mxr]))-dp),-dp)
    x = np.zeros(p.size)
    y = np.zeros(p.size)

    for mix0 in mixrat:
        for k in range(p.size):
            mix = Td(mix0/1000.,p[k]/1000.)-T0 
            y[k] = skewty(p[k])
            xtmp = skewtx(mix,y[k])
            if(xtmp >= x00 and xtmp <= x11):
                x[k] = xtmp
            else:
                x[k] = -9999

        doPlot = np.where(x != -9999)
        if(len(doPlot[0]) > 0):
            pl.plot(x[doPlot],y[doPlot],color=c5,dashes=[3,2])
            pl.text(x[doPlot][-1],y[doPlot][-1]+vs,int(mix0),color=c5,ha='center',va='bottom',backgroundcolor='w')

    """
    6. Add sounding data
    """
    # 6.1 Temperature and dewpoint temperature
    if(si.T.size > 0):
        y = skewty(si.p)
        x1 = skewtx(si.T-T0,y)
        x2 = skewtx(si.Td-T0,y)

        pl.plot(x1,y,'-',color=cT,linewidth=2)
        pl.plot(x2,y,'-',color=cTd,linewidth=2)

    # 6.2 Add height labels to axis
    if(si.z.size > 0 and si.z.size==si.p.size):
        y = skewty(si.p)
        p_last = 1e9 
        for k in range(si.z.size):
            if(y[k] <= y11 and np.abs(si.p[k]-p_last) > dp_label):
                pl.text(x11+hs,y[k],str(int(si.z[k]))+'m',color=c2,ha='right',va='center',size=7,backgroundcolor='w')
                p_last = si.p[k]

    # 6.2 Wind barbs
    if(si.u.size > 0):
        y   = skewty(si.p)
        u   = si.u * 1.95
        v   = si.v * 1.95 
        xb  = x11 + 9*hs 

        p_last = 1e9 
        for k in range(si.z.size):
            if(y[k] <= y11 and np.abs(si.p[k]-p_last) > dp_label):
                pl.barbs(xb,y[k],u[k],v[k],length=5.2,linewidth=0.5,pivot='middle')  
                p_last = si.p[k]

    """
    7. Lauch parcel
    """
    if(si.parcel == True):

        # HACK BvS
        si.ps = si.p[0]

        # Iteratively find Lifting Condensation Level (LCL)
        zlcl  = 1000. 
        RHlcl = 0.9
        niter = 0;
        while(((RHlcl <= 0.9999) or (RHlcl >= 1.0001)) and niter < 30): 
            zlcl   += (1.-RHlcl) * 1000.;
            Plcl    = si.ps / np.exp((g * zlcl)/(Rd * si.Ts));
            Tlcl    = si.Ts * exner(Plcl,si.ps) # pow(Ps / Plcl,Rd / cp);
            esatlcl = 0.611e3 * np.exp((Lv / Rv) * ((1. / 273.15)-(1. / Tlcl)));
            elcl    = si.rs * Plcl / 0.622;
            RHlcl   = elcl / esatlcl;
            niter  += 1;

        # Surface values        
        psurfs = skewty(si.ps)
        Tsurfs = skewtx(si.Ts-T0,psurfs)
        Tdsurfs = skewtx(Td(si.rs,si.ps/1000.)-T0,psurfs)
        pl.scatter(Tsurfs,psurfs,facecolor='none')
        pl.scatter(Tdsurfs,psurfs,facecolor='none')

        # Dry adiabat + mixing ratio up to LCL
        plcls  = skewty(Plcl)
        Tlcls  = skewtx(Tlcl-T0,plcls)
        pl.scatter(Tlcls,plcls,facecolor='none')
        pl.plot([Tsurfs,Tlcls],[psurfs,plcls],'k-') 
        pl.plot([Tdsurfs,Tlcls],[psurfs,plcls],'k-') 
        pl.plot([Tlcls-5*hs,Tlcls+5*hs],[plcls,plcls],'k',dashes=[2,1])
        pl.text(Tlcls+6*hs,plcls,'LCL',ha='left',va='center',size=8)

        # Moist adiabat from LCL upwards
        Ths    = si.Ts / exner(si.ps) 
        thw0   = dsatlftskewt(Ths-T0,Plcl)+T0 
        thwlcl = dsatlftskewt(thw0-T0,Plcl)+T0
        thw0  -= thwlcl - Tlcl

        p = np.arange(Plcl,ptop,-dp)
        x = np.zeros(p.size)
        y = np.zeros(p.size)
        for k in range(p.size):
            thw = dsatlftskewt(thw0-T0,p[k])
            y[k] = skewty(p[k])
            x[k] = skewtx(thw,y[k])

        pl.plot(x,y,'k-')


    """
    6. Finish diagram
    """ 
    pl.xticks([])
    pl.yticks([])
    pl.box('off')
    pl.xlim(x00-0.1,x11+16*hs)
    pl.ylim(y00-0.1,y11+0.1)

    # draw axis
    pl.plot([x00,x11],[y00,y00],'k-')
    pl.plot([x00,x00],[y00,y11],'k-')

    pl.figtext(0.5,0.05,'Temperature (C)',ha='center')
    pl.figtext(0.015,0.52,'Pressure (hPa)',rotation=90)
    label = 'Skew-T log-P, %s, %s'%(si.name,si.time) 
    pl.figtext(0.5,0.97,label,ha='center')

    img = image.imread('data/olga_lr.png')
    w=650;h=600
    pl.figimage(img,10,6)
    pl.figtext(0.08,0.013,'Open Limited-area Gliding Analysis. 6 x 6 km GFS-initiated WRF-ARW forecast [olga.vanstratum.com]',size=7,ha='left')

"""
Just for testing..
"""
if __name__ == "__main__":
    pl.close('all')
    sset = skewt_input()
    skewtlogp(sset)
    sset.stype = 1
    skewtlogp(sset)

