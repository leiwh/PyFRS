#==================================Astrophysics Fundaments=================================
#Credit: Weihua Lei
#Institue for Particle Physics and Astrophysics, Huazhong University of Science and Technology (HUST), Wuhan, China
#Email: leiwh@hust.edu.cn
#Homepage: http://www.physics.unlv.edu/~leiwh/pygrb/
#
#  Description:
#Created on Tue Aug 16 19:39:07 2011, University of Nevada Las Vegas (UNLV), USA
#Current version: 2014-11-25
#
#
#==============================================================================

import math
import cgs
import numpy as np
from scipy.integrate import quad
from cosmocalc import cosmocalc

# Hubble constant, cm/s/cm
H_0=72.0 *1.0e5/(1.0e6 *cgs.pc); 
def flumdist(x):
   return 1.0/math.sqrt(0.3 *(1.0+x)**3.0 + 0.7)

def lumdist(z=None):
    if (z is None):
        ss="Syntax: lumdist(z)"
    else:    
        finte=quad(flumdist,0,z)[0]
        ss=cgs.c*(1.0+z)/H_0 *finte
    return ss

def LUD(z):
    return cosmocalc(z)['DL_cm']


#Klein-Nishina scattering cross section
def sigmaKN(x):
   if (x > 1.e-4):
      s1=(1.+x)/x**3. *(2.*x*(1.+x)/(1.+2.*x)-np.log(1.+2.*x) )
      s2=1./(2.*x) *np.log(1+2.*x) - (1.+3.*x)/(1.+2.*x)**2.
      s=0.75*cgs.sigmaT*(s1+s2)
   else:
      s=cgs.sigmaT
   return s


def band(UBV,key):
    if (key is None):
        ss= "Syntax: band(UBV,key); About: Standard apparent magnitudes and fluxes; flux in Jy, lambda in um; nu in Hz"
    else:
        flux={
            'U': [0.36, 1810],
            'B': [0.22, 4260],
            'V': [0.55, 3640],
            'R': [0.64, 3080],
            'I': [0.79, 2550],
            'J': [1.26, 1600],
            'H': [1.6, 1080],
            'K': [2.22, 670],
            'g': [0.52, 3730],
            'r': [0.67, 4490],
            'i': [0.79, 4760],
            'z': [0.91, 4810]
        }
        ikey = {'lambda': 0, 'flux': 1}
        if (key == 'nu'):
            ss= 1.e4*cgs.c/flux[UBV][ikey['lambda']]
        else:
            ss = flux[UBV][ikey[key]]
    return ss

def mag(UBV, m=None):
    if (m is None):
        ss= "Syntax: mag(UBV,m); About: aparent magnitudescy, unit: erg cm^-2 Hz^-1 s^-1"
    else:
        ss=band(UBV,'flux')*cgs.Jy*100.0**(-m/5.)
    return ss


def m2MAG(d, m=None):
    if (m is None):
        ss= "Syntax: m2MAG(d, m); apparent mag (m) to absoulte mag (M); About: distance d in pc"
    else:
        ss= m - 5.* np.log10(d/10.)
    return ss


def Msmbh(Lbul=None):
    if (Lbul is None):
        ss= "Syntax: MBH(Lbulge); About: BH mass of galatic center SMBH with M-Lbul relation, unit: Msun; for Lbul: Lsun"
    else:
        ss=10.**(8.07+1.26*(np.log10(Lbul)-10.) )
    return ss

def darcsec(z=None):
    if (z is None):
        ss= "Syntax: darcsec(z); About: distance for 1 arcsec for object with known redshift, where z=0 for objects in Milky Way, unit: kpc"
    else:
        if z==0.0:
            dl = cgs.Dsun
        else:
            dl= LUD(z)
            da=dl/(1.+z)**2.
        ss= da*cgs.pi/180./3600./cgs.kpc 
    return ss

    
def nu_cy(B=None):
    if (B is None):
        ss= "Syntax: nu_cy(B); About: cyclotron frequency"
    else:
        ss=cgs.eq*B/(2.*cgs.pi*cgs.me*cgs.c)
    return ss
    
def r_A(B=None, Mns=1.4, rns=1e6,mdot=5e-17):
    if (B is None):
        ss= "Syntax: r_A(B,Mns, rns,mot); About: Alfven radius"
    else:
        ss=(B**4.*rns**12./(2.*cgs.G*Mns*cgs.Msun*3.*mdot**2.))**(1/7.)
    return ss

def t_ff(rho=None):
    if (rho is None):
        ss= "Syntax: t_ff(rho); About: free-fall time for a cloud that has satisfied the Jeans criterion"
    else:
        ss=(3.*cgs.pi/(32.*cgs.G*rho))**0.5
    return ss
    
def R_J(rho=None, T=50.):
    mu=1.
    if (rho is None):
        ss= "Syntax: R_J(rho, T); About: Jeans length"
    else:
        ss=(15.*cgs.k *T/(4.*cgs.pi*cgs.G*mu*cgs.mH*rho))**0.5
    return ss    


def gamma(beta=None):
    if (beta is None):
        ss= "Syntax: gamma(beta); About: velocity to Gamma"
    else:
        ss=(1/(1-beta**2))**0.5
    return ss    

def beta(gamma=None):
    if (beta is None):
        ss= "Syntax: beta(gamma); About: gamma to velocity"
    else:
        ss=(1-1/gamma**2)**0.5
    return ss    
