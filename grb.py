#==================================PyGRB Toolkit===--==============================
#  A Python tool for GRB model Fit
#Credit: Weihua Lei
#Institue for Particle Physics and Astrophysics, Huazhong University of Science and Technology (HUST), Wuhan, China
#Email: leiwh@hust.edu.cn
#Homepage: http://www.physics.unlv.edu/~leiwh/pygrb/
#
#  Description:
#Can be applied to GRB, TDE, etc.
#Models: forward shock, energy injection, two-component jet, SSC, IC by X-ray photons from internal shock (XIC)
#Created on Tue Aug 16 19:39:07 2011, University of Nevada Las Vegas (UNLV), USA
#Current version: 2014-11-25
#
#  Reference:
#-------afterglow dynamics and emission:
#1. A complete reference of the analytical synchrotron external shock models of gamma-ray bursts, He Gao, Wei-Hua Lei, Yuan-Chuan Zou, Xue-Feng Wu, Bing Zhang,2013, New Astronomy Review, 57, 141-190
#2. Overall Evolution of Jetted Gamma-Ray Burst Ejecta, Y. F. Huang, L. J. Gou, Z. G. Dai, T. Lu, 2000, ApJ, 543, 90-96
#-------central engine models:
#1. Hyperaccreting Black Hole as Gamma-Ray Burst Central Engine. I. Baryon Loading in Gamma-Ray Burst Jets, Wei-Hua Lei, Bing Zhang, En-Wei Liang, 2013, ApJ, 765, 125
#2. Giant X-Ray Bump in GRB 121027A: Evidence For Fall-Back Disk Accretion, Xue-Feng Wu, Shu-Jin Hou, Wei-Hua Lei, 2013, ApJ, 765, 36 
#
#
#  These papers used our code:
#QUASI-PERIODIC VARIATIONS IN X-RAY EMISSION AND LONG-TERM RADIO OBSERVATIONS: EVIDENCE FOR A TWO-COMPONENT JET IN Sw J1644+57, Jiu-Zhou Wang, Wei-Hua Lei, Ding-Xiong Wang, Yuan-Chuan Zou, Bing Zhang, He Gao, Chang-Yin Huang, 2014, ApJ, 788, 32
#
#==============================================================================

import math
import cgs,astro
import numpy as np
from scipy.special import kv
from scipy.integrate import quad

#Doppler ratio for off-beaming
#off-beaming angle
def theta_off(GM,thetaobs,thetaj):
    thetacone=1./GM
    if (thetaj>=thetacone):
        thetajet=thetaj
    else:
        thetajet=thetacone

    if (thetaobs <= thetajet):
        theta_view=0.
    else:
        theta_view=(thetaobs - thetajet)
    
    return theta_view


def aoff(GM,theta_view):
    beta=math.sqrt(1.-1./GM**2.0)
    return (1.-beta)/(1.-beta*math.cos(theta_view) )



#Off-beaming correction factor to Fv
def FvOff(GM,thetaobs,thetaj):
    thetacone=1./GM
    theta_view=theta_off(GM,thetaobs,thetaj)
    fview=aoff(GM,theta_view)    
    if (thetaj<=thetacone):
        if ((thetaobs>=0.) and (thetaobs<=-thetaj+thetacone)):
            fviewF=1.
        elif ((thetaobs>-thetaj+thetacone) and (thetaobs<=thetacone)):
            fviewF=1.-(thetaobs+thetaj-thetacone)/2./thetaj
        elif ((thetaobs>thetacone) and (thetaobs<2.*thetacone)):
            fviewF=fview**2. *fview*(1.+(GM* theta_view)**2.)/2.
        elif (thetaobs>=2.*thetacone):
            fviewF=fview**3.
    else:
        if ((thetaobs>=0.) and (thetaobs<=thetaj-thetacone)):
            fviewF=1.
        elif ((thetaobs>thetaj-thetacone) and (thetaobs<=thetaj)):
            fviewF=1.-GM*(thetaobs-thetaj+thetacone)/2.
        elif ((thetaobs>thetaj) and (thetaobs<2.*thetaj)):
            fviewF=fview**2. *fview*(1.+(GM* theta_view)**2.)/2.
        elif (thetaobs>=2.*thetaj):
            fviewF=fview**3. *(1.+(GM* thetaj)**2.)/2.
    return fviewF




#time for nu_m cross
def tm(eb,ee,E52,nu15):
    return 0.69*eb**(1/3.)*ee**(4/3.)*E52**(1/3.)*nu15**(-2/3.)

#time for nu_c cross    
def tc(eb,E52,n,nu15):
    return 7.3e-6*eb**(-3.)*n**(-2.)*E52**(-1.)*nu15**(-2.)

#band function
def NE(E,alpha,beta,E0):
    if (E<(alpha-beta)*E0):
        s=E**alpha *math.exp(-E/E0)
    else:
        s=((alpha-beta)*E0)**(alpha-beta)*E**beta*math.exp(beta-alpha)
    return s

#bulk Lorentz factor Gamma as a function of R
#valid only for deceraction phase, no energy injection, ISM case 
def Gammat2(Gamma0,R,R0):
    return Gamma0*(R/R0)**(-3/2.)
    
#bulk Lorentz factor Gamma as a function of R
#valid for both pre and after deceraction radius, no energy injection, ISM case 
def Gammat1(Gm0,n, Eiso, R,R0):
    eps=1.e-6
    A=4.*cgs.pi*R**3.*n*cgs.mp/3.
    B=Eiso/(Gm0*cgs.c**2.)-A
    C=Eiso/cgs.c**2.
    ABC=4.*A*C/B**2.
    if (ABC <eps):
        Gmt=Gm0
    else:
        Gmt=(-B+math.sqrt(B**2.* (1.+ABC)))/(2.*A)
    return Gmt

#bulk Lorentz factor Gamma as a function of R
#valid for both pre and after deceraction radius, no energy injection, ISM case
def Gammat(Gm0,n, Eiso, R,R0):
    eps=1.e-6
    A=4.*cgs.pi*R**3.*n*cgs.mp/3.
    B=Eiso/(Gm0*cgs.c**2.)
    C=Eiso/cgs.c**2.+A
    ABC=4.*A*C/B**2.
    if (ABC <eps):
        Gmt=Gm0+A/B
    else:
        Gmt=(-B+math.sqrt(B**2.* (1.+ABC)))/(2.*A)
    return Gmt
 
#Deceration radius   
def Rdec(Gm0,n, Eiso):
    return (3.*Eiso/(Gm0**2.* cgs.c**2. *2.*cgs.pi*n*cgs.mp))**(1/3.)

#Deceration time
def Tdec(Gm0,n,Eiso):
    return Rdec(Gm0,n,Eiso)/2./Gm0**2/cgs.c

#Sedov radius
def Rsed(n, Eiso):
    return (3.*Eiso/(cgs.c**2 *4.*cgs.pi*n*cgs.mp))**(1/3.)

#Comoving B, version 1: with input parameters: epsilon_B, Gamma, n 
def Bco(epsilon_B, Gamma0, n1):
    return math.sqrt(8.*cgs.pi*(4.*Gamma0+3.)*(Gamma0-1.)*cgs.mp*epsilon_B*n1)*cgs.c

#Comoving B, version 2: with input parameters: epsilon_B, e
def Bco2(epsilon_B, ei):
    return math.sqrt(8.*cgs.pi* epsilon_B* ei)

def n2_c(GM21,n1):
    gmcat = (4.*GM21+1.)/(3.*GM21)
    n2=(gmcat*GM21+1.)/(gmcat-1.)*n1
    return n2


def e2(GM21,n1):
    gmcat = (4.*GM21+1.)/(3.*GM21)
    if (GM21-1.) > 1.e-3:
#--------For Relativistic case        
#        gmcat=4./3.
#        n2n1=(gmcat*GM21+1.)/(gmcat-1.)
        GMm1 = GM21 -1.
    else:
#--------Newtonian case        
#        gmcat=5./3.
        beta=math.sqrt(1.-1./GM21**2.0)
        GMm1 = 0.5*(beta)**2
    n2n1=(gmcat*GM21+1.)/(gmcat-1.)
    ei= n2n1* GMm1 *n1 *cgs.mp *cgs.c**2
    return ei

#--------------------------------------------------------------------------
#---------------------------------gamma--------------------------------------
#--------------------------------------------------------------------------

#max electron Lorentz factor: gamma_m
def gamma_Max(Bc):
    zeta = 1.
    return math.sqrt(6.*cgs.pi*cgs.qe/(cgs.sigmaT* Bc* zeta))

#minimum electron Lorentz factor: gamma_m
def gamma_m(epsilon_e,GM,p):   
#---for p>2
    g_p=(p-2.)/(p-1.)
    npne=1.0
    gamma_ave= epsilon_e*  (GM-1.)*npne* cgs.mp/cgs.me  
    return g_p* gamma_ave

def gamma_m2(epsilon_e,GM,p, g_M):
#---for 1<p<2, p!=1, p!=2, include Newtonian case
    if p > 2.:
        if (GM-1.) > 1.e-3:
            gsm = gamma_m(epsilon_e,GM,p)
        else:
            beta=math.sqrt(1.-1./GM**2.0)
            g_p=(p - 2.)/(p - 1.)
            npne=1.0
            gamma_ave=epsilon_e* 0.5*(beta)**2 *npne* cgs.mp/cgs.me 
            gsm = g_p* gamma_ave
    elif 1.< p < 2.:
        pp = (2.-p)/(p-1.)
        npne = 1.
        mpme = cgs.mp/cgs.me
        GMm1 = GM -1.
        if (GM-1.) <= 1.e-3:
            beta=math.sqrt(1.-1./GM**2.0)
            GMm1 = 0.5*(beta)**2
        gsm = (epsilon_e *GMm1 *pp*mpme*npne/g_M**(2.-p) )**(1./(p-1.))
    return gsm

    
#cooling electron Lorentz factor: gamma_c    
def gamma_c(Bc,GM,t,z):
    return 6.*cgs.pi*cgs.me*cgs.c/(cgs.sigmaT* Bc**2* GM* t/(1+z))

#get nu in lab frame by given electron Lorentz factor 
def nu_gme(GM, gm_e, Bc, z):
    return GM* gm_e**2 *cgs.qe *Bc/(2.*cgs.pi *cgs.me *cgs.c)/(1.+z)
#    return GM* gm_e**2 *cgs.qe *Bc*3./(4.*cgs.pi *cgs.me *cgs.c)/(1.+z)

#get maximum frequence in Lab frame: nu_Max
def nu_Max(GM, z):
    return 3.* GM *cgs.qe**2/(2.*cgs.pi *cgs.me *cgs.c)/(1.+z)/(cgs.sigmaT)

def Fvmax_analy(epsilon_B, n, E52, D28):
    return 1.1e5*epsilon_B**(1/2.)*E52*n**(1/2.)*D28**(-2.)

def Pvmax(GM, Bc):
    return cgs.me*cgs.c**2.0*cgs.sigmaT*GM*Bc/(3.*cgs.qe)

def Fvmax(Pvm,n1,R,D28,z):
    Ne=4.*cgs.pi*R**3.*n1/3
    return (1+z)*Ne*Pvm/(4.*cgs.pi*(D28*1e28)**2.)
    
def Fvmax_IC(Fvm, n1, R):
    x0=0.4714 #sqrt(2)/3
    return cgs.sigmaT *n1 *R *Fvm* x0

def Fvmax_IC0(Fvm, n1, R):
    return cgs.sigmaT *n1 *R *Fvm/3.

def Fvmax3(m,epsilon_B,Gamma0, n1,D28):
    Ne=m/cgs.mp
    B=Bco(epsilon_B,Gamma0,n1)
    Pvm=Pvmax(Gamma0,B)
    return Ne*Pvm/(4.*cgs.pi*(D28*1e28)**2.)/cgs.uJy


#--------------------------------------------------------------------------
#---------------------------------nu--------------------------------------
#--------------------------------------------------------------------------
#time evolution of nu_c and nu_m
#valid only for deceration phase 
def nu_c(epsilon_B, n, E52, td):
    return 2.7e12* epsilon_B**(-3/2.)*E52**(-1/2.)*td**(-1/2.)/n

def nu_m(epsilon_B, epsilon_e, E52, td):
    return 5.7e14* epsilon_B**(1/2.)*epsilon_e**2.*E52**(1/2.)*td**(-3/2.)




#nu=GM* gm_e**2 *cgs.qe *Bc/(2.*cgs.pi *cgs.me *cgs.c)/(1.+z)
#the first version of nu_a
def nu_a1(epsilon_B,p,n1,vm,vc,Gamma0,R,z):
    B=Bco(epsilon_B,Gamma0,n1)
    gamma_nu=(2.*cgs.pi*cgs.me*cgs.c*(1.+z)/(Gamma0*cgs.qe*B))**0.5
    gamma_m= vm**0.5 *gamma_nu
    gamma_c= vc**0.5 *gamma_nu
    if (vm<vc):
        Cg1=4.*Gamma0*n1*(p-1.)*gamma_m**(p-1.)
        alpha=1.0e4*(8.4e6)**(p/2.)*Cg1*B**(p+2.)/2.*math.gamma((3.*p+2.)/12.)*math.gamma((3.*p+22.)/12.)*R/Gamma0
        tau_vm=alpha*(vm*(1.+z)/Gamma0)**(-(p+4.)/2.)
        tau_vc=alpha*(vc*(1.+z)/Gamma0)**(-(p+4.)/2.)
        if (tau_vm <1.):
            Cg=4.*Gamma0*n1*(p-1.)*gamma_m**(p-1.)
            va=(136.*(p+2.)/(p+2./3.)*Cg*B**(2./3.)*gamma_m**(-(p+2./3.)) *R/Gamma0)**(3./5.) *Gamma0/(1.+z)
        else:
            if (tau_vc <1.):
                Cg=4.*Gamma0*n1*(p-1.)*gamma_m**(p-1.)
#                va= alpha**(2./(p+4.))*Gamma0/(1.+z)
                pj=p
                va=(1.0e4*(8.4e6)**(pj/2.)*Cg*B**((pj+2.)/2.)*math.gamma((3.*pj+2.)/12.)*math.gamma((3.*pj+22.)/12.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)
            else:
                Cg=4.*Gamma0*n1*(p-1.)*gamma_m**(p-1.)*gamma_c
                pj=p+1.
                va=(1.0e4*(8.4e6)**(pj/2.)*Cg*B**((pj+2.)/2.)*math.gamma((3.*pj+2.)/12.)*math.gamma((3.*pj+22.)/12.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)
    else:
        Cg1=4.*Gamma0*n1*gamma_c
        alpha=1.0e4*8.4e6*Cg1*B**2.*math.gamma((3.*2.+2.)/12.)*math.gamma((3.*2.+22.)/12.)*R/Gamma0
        tau_vm=alpha*(vm*(1.+z)/Gamma0)**(-3.)
        tau_vc=alpha*(vc*(1.+z)/Gamma0)**(-3.)
        if (tau_vc <1.):
            Cg=4.*Gamma0*n1*gamma_c
            va=(136.*(2.+2.)/(2.+2./3.)*Cg*B**(2./3.)*gamma_c**(-(2.+2./3.))*R/Gamma0)**(3/5.)*Gamma0/(1.+z)
        else:
            if (tau_vm <1.):
                Cg=4.*Gamma0*n1*gamma_c
#                va= alpha**(1/3.)*Gamma0/(1.+z)
                pj=2.
                va=(1.0e4*(8.4e6)**(pj/2.)*Cg*B**((pj+2.)/2.)*math.gamma((3.*pj+2.)/12.)*math.gamma((3.*pj+22.)/12.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)    
            else:
                Cg=4.*Gamma0*n1*gamma_m**(p-1.)*gamma_c
                pj=p+1.
                va=(1.0e4*(8.4e6)**(pj/2.)*Cg*B**((pj+2.)/2.)*math.gamma((3.*pj+2.)/12.)*math.gamma((3.*pj+22.)/12.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)    
    return va
    
#the second version of nu_a
def nu_a2(epsilon_B,p,n1,vm,vc,Gamma0,R,z):
    B=Bco(epsilon_B,Gamma0,n1)
    gamma_nu=(2.*cgs.pi*cgs.me*cgs.c*(1.+z)/(Gamma0*cgs.qe*B))**0.5
    gamma_m= vm**0.5 *gamma_nu
    gamma_c= vc**0.5 *gamma_nu
    if (vm<vc):
        Cg=4.*Gamma0*n1*(p-1.)*gamma_m**(p-1.)
        va1=(136.*(p+2.)/(p+2./3.)*Cg*B**(2./3.)*gamma_m**(-(p+2./3.)) *R/Gamma0)**(3./5.) *Gamma0/(1.+z)
        Cg=4.*Gamma0*n1*(p-1.)*gamma_m**(p-1.)
        pj=p
        va2=(1.0e4*(8.4e6)**(pj/2.)*Cg*B**((pj+2.)/2.)*math.gamma((3.*pj+2.)/12.)*math.gamma((3.*pj+22.)/12.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)
        Cg=4.*Gamma0*n1*(p-1.)*gamma_m**(p-1.)*gamma_c
        pj=p+1.
        va3=(1.0e4*(8.4e6)**(pj/2.)*Cg*B**((pj+2.)/2.)*math.gamma((3.*pj+2.)/12.)*math.gamma((3.*pj+22.)/12.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)
        if (va1 <= va2):
            va=va1
        else:
            if (va2 <= va3):
                va=va2
            else:
                va=va3
    else:
        Cg=4.*Gamma0*n1*gamma_c
        va1=(136.*(2.+2.)/(2.+2./3.)*Cg*B**(2./3.)*gamma_c**(-(2.+2./3.))*R/Gamma0)**(3/5.)*Gamma0/(1.+z)
        Cg=4.*Gamma0*n1*gamma_c
        pj=2.
        va2=(1.0e4*(8.4e6)**(pj/2.)*Cg*B**((pj+2.)/2.)*math.gamma((3.*pj+2.)/12.)*math.gamma((3.*pj+22.)/12.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)    
        Cg=4.*Gamma0*n1*gamma_m**(p-1.)*gamma_c
        pj=p+1.
        va3=(1.0e4*(8.4e6)**(pj/2.)*Cg*B**((pj+2.)/2.)*math.gamma((3.*pj+2.)/12.)*math.gamma((3.*pj+22.)/12.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)    
        if (va1 <= va2):
            va=va1
        else:
            if (va2 <= va3):
                va=va2
            else:
                va=va3

    return va

#the third version of nu_a
#nu_a is not smoothly connected    
def nu_a3(epsilon_B,p,n1,vm,vc,Gamma0,R,z):
    B=Bco(epsilon_B,Gamma0,n1)
    gamma_nu=(2.*cgs.pi*cgs.me*cgs.c*(1.+z)/(Gamma0*cgs.qe*B))**0.5
    gamma_m= vm**0.5 *gamma_nu
    gamma_c= vc**0.5 *gamma_nu
    if (vm<vc):
        Cg=4.*Gamma0*n1*(p-1.)*gamma_m**(p-1.)
        va1=(136.*(p+2.)/(p+2./3.)*Cg*B**(2./3.)*gamma_m**(-(p+2./3.)) *R/Gamma0)**(3./5.) *Gamma0/(1.+z)
        Cg=4.*Gamma0*n1*(p-1.)*gamma_m**(p-1.)
        pj=p
        va2=(1.0e4*(8.4e6)**(pj/2.)*Cg*B**((pj+2.)/2.)*math.gamma((3.*pj+2.)/12.)*math.gamma((3.*pj+22.)/12.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)
        Cg=4.*Gamma0*n1*(p-1.)*gamma_m**(p-1.)*gamma_c
        pj=p+1.
        va3=(1.0e4*(8.4e6)**(pj/2.)*Cg*B**((pj+2.)/2.)*math.gamma((3.*pj+2.)/12.)*math.gamma((3.*pj+22.)/12.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)
        if (va1 <= vm):
            va=va1
        else:
            if (va2 <= vc):
                va=va2
            else:
                va=va3
    else:
        Cg=4.*Gamma0*n1*gamma_c
        va1=(136.*(2.+2.)/(2.+2./3.)*Cg*B**(2./3.)*gamma_c**(-(2.+2./3.))*R/Gamma0)**(3/5.)*Gamma0/(1.+z)
        Cg=4.*Gamma0*n1*gamma_c
        pj=2.
        va2=(1.0e4*(8.4e6)**(pj/2.)*Cg*B**((pj+2.)/2.)*math.gamma((3.*pj+2.)/12.)*math.gamma((3.*pj+22.)/12.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)    
        Cg=4.*Gamma0*n1*gamma_m**(p-1.)*gamma_c
        pj=p+1.
        va3=(1.0e4*(8.4e6)**(pj/2.)*Cg*B**((pj+2.)/2.)*math.gamma((3.*pj+2.)/12.)*math.gamma((3.*pj+22.)/12.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)    
        if (va1 <= vc):
            va=va1
        else:
            if (va2 <= vm):
                va=va2
            else:
                va=va3
    return va

#the fourth of nu_a
#modified for smooth connection between transition
def nu_a(epsilon_B,p,n1,vm,vc,Gamma0,R,z):
    B=Bco(epsilon_B,Gamma0,n1)
    c0=(2.**(4/3.)*math.gamma(1/3.))**(-1.) * 4.*cgs.pi**2. *cgs.qe
    cnugm=2.*cgs.pi*cgs.me*cgs.c/cgs.qe
    gamma_nu=(cnugm*(1.+z)/B/Gamma0)**0.5
    gamma_m= vm**0.5 *gamma_nu
    gamma_c= vc**0.5 *gamma_nu
    if (vm<vc):
        Cg=4.*Gamma0*n1*(p-1.)*gamma_m**(p-1.)
        pj=p
        va1=((c0/B)*(pj+2.)/(pj+2./3.)*Cg*(cnugm/B)**(-5./3.)*gamma_m**(-(pj+2./3.)) *R/Gamma0)**(3./5.) *Gamma0/(1.+z)
        Cg=4.*Gamma0*n1*(p-1.)*gamma_m**(p-1.)
        pj=p
        va2=((c0/B)*(pj+2.)/(pj+2./3.)*Cg*(cnugm/B)**(-(pj+4.)/2.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)
        Cg=4.*Gamma0*n1*(p-1.)*gamma_m**(p-1.)*gamma_c
        pj=p+1.
        va3=((c0/B)*(pj+2.)/(pj+2./3.)*Cg*(cnugm/B)**(-(pj+4.)/2.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)
        if (va1 <= vm):
            va=va1
        else:
            if (va2 <= vc):
                va=va2
            else:
                va=va3
    else:
        Cg=4.*Gamma0*n1*gamma_c
        pj=2.
        va1=((c0/B)*(pj+2.)/(pj+2./3.)*Cg*(cnugm/B)**(-5./3.)*gamma_c**(-(pj+2./3.))*R/Gamma0)**(3/5.)*Gamma0/(1.+z)
        Cg=4.*Gamma0*n1*gamma_c
        pj=2.
        va2=((c0/B)*(pj+2.)/(pj+2./3.)*Cg*(cnugm/B)**(-(pj+4.)/2.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)    
        Cg=4.*Gamma0*n1*gamma_m**(p-1.)*gamma_c
        pj=p+1.
        va3=((c0/B)*(pj+2.)/(pj+2./3.)*Cg*(cnugm/B)**(-(pj+4.)/2.)*R/Gamma0)**(2./(pj+4.))*Gamma0/(1.+z)    
        if (va1 <= vc):
            va=va1
        else:
            if (va2 <= vm):
                va=va2
            else:
                va=va3
    return va

#the fivth and current version of nu_a with given B
#modified for smooth connection between transition
def nu_aB(Bc,p,n1,vm,vc,Gamma,R,z):
    B=Bc
#    c0=(2.**(4/3.)*math.gamma(1/3.))**(-1.) * 4.*cgs.pi**2. *cgs.qe
    c0= 2.8071453296870878e-09
#    cnugm=2.*cgs.pi*cgs.me*cgs.c/cgs.qe
    cnugm=3.575013703788261e-07
    gamma_nu=math.sqrt(cnugm*(1.+z)/B/Gamma)
    gamma_m= math.sqrt(vm) *gamma_nu
    gamma_c= math.sqrt(vc) *gamma_nu
    if (vm<vc):
        Cg=4.*Gamma*n1*(p-1.)*gamma_m**(p-1.)
        pj=p
        va1=((c0/B)*(pj+2.)/(pj+2./3.)*Cg*(cnugm/B)**(-5./3.)*gamma_m**(-(pj+2./3.)) *R/Gamma)**(3./5.) *Gamma/(1.+z)
#        Cg=4.*Gamma*n1*(p-1.)*gamma_m**(p-1.)
#        pj=p
        va2=((c0/B)*(pj+2.)/(pj+2./3.)*Cg*(cnugm/B)**(-(pj+4.)/2.)*R/Gamma)**(2./(pj+4.))*Gamma/(1.+z)
        Cg=4.*Gamma*n1*(p-1.)*gamma_m**(p-1.)*gamma_c
        pj=p+1.
        va3=((c0/B)*(pj+2.)/(pj+2./3.)*Cg*(cnugm/B)**(-(pj+4.)/2.)*R/Gamma)**(2./(pj+4.))*Gamma/(1.+z)
        if (va1 <= vm):
            va=va1
        else:
            if (va2 <= vc):
                va=va2
            else:
                va=va3
    else:
        Cg=4.*Gamma*n1*gamma_c
        pj=2.
        va1=((c0/B)*(pj+2.)/(pj+2./3.)*Cg*(cnugm/B)**(-5./3.)*gamma_c**(-(pj+2./3.))*R/Gamma)**(3/5.)*Gamma/(1.+z)
#        Cg=4.*Gamma*n1*gamma_c
#        pj=2.
        va2=((c0/B)*(pj+2.)/(pj+2./3.)*Cg*(cnugm/B)**(-(pj+4.)/2.)*R/Gamma)**(2./(pj+4.))*Gamma/(1.+z)    
        Cg=4.*Gamma*n1*gamma_m**(p-1.)*gamma_c
        pj=p+1.
        va3=((c0/B)*(pj+2.)/(pj+2./3.)*Cg*(cnugm/B)**(-(pj+4.)/2.)*R/Gamma)**(2./(pj+4.))*Gamma/(1.+z)    
        if (va1 <= vc):
            va=va1
        else:
            if (va2 <= vm):
                va=va2
            else:
                va=va3
    return va

#--------------------------------------------------------------------------
#---------------------------------Syn--------------------------------------

#-------------------------Trail version of Fv---------------------
#no nu_a
def Fv_trail(v, vm, vc, vmax, Fvm, p):
    if (vc<=vm):
#------fast cooling-------
        if (v<=vc):
            Fnu=(v/vc)**(1/3.)* Fvm
        elif (v>vc) and (v<=vm):
            Fnu=(v/vc)**(-1/2.)*Fvm
        elif (v>vm) and (v<=vmax):
            Fnu=(vm/vc)**(-1/2.)*(v/vm)**(-p/2.)*Fvm
        elif (v>vmax):
            Fnu=(vm/vc)**(-1/2.)*(vmax/vm)**(-p/2.)*Fvm * np.exp(1.-v/vmax)
    else:
#------slow cooling------
        if (v<=vm):
            Fnu=(v/vm)**(1/3.)* Fvm
        elif (v>vm) and (v<=vc):
            Fnu=(v/vm)**(-(p-1.)/2.)*Fvm
        elif (v>vc) and (v<=vmax):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
        elif (v>vmax):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(vmax/vc)**(-p/2.)*Fvm * np.exp(1.-v/vmax)
    return Fnu


#------------------------Simple braken power law---------------------------
#Fv in fast cooling case
def Fvfast_sbpl(v,va,vm,vc,Fvm,p):
    if (va<=vc):
        if (v<=va):
            Fnu=(va/vc)**(1/3.) *(v/va)**2. *Fvm
        elif (v>va) and (v<=vc):
            Fnu=(v/vc)**(1/3.)* Fvm
        elif (v>vc) and (v<=vm):
            Fnu=(v/vc)**(-1/2.)*Fvm
        elif (v>vm):
            Fnu=(vm/vc)**(-1/2.)*(v/vm)**(-p/2.)*Fvm
    elif (va>vc) and (va<=vm):
        if (v<=vc):
            Fnu=(vc/va)**3. *(v/vc)**2. *Fvm
        elif (v>vc) and (v<=va):
            Fnu=(va/vc)**(-1/2.) *(v/va)**2.5 *Fvm
        elif (v>va) and (v<=vm):
            Fnu=(v/vc)**(-1/2.)*Fvm
        elif (v>vm):
            Fnu=(vm/vc)**(-1/2.)*(v/vm)**(-p/2.)*Fvm    
    elif (va>vm):
        if (v<=vc):
            Fnu=(vc/va)**3. *(vm/va)**((p-1.)/2.) *(v/vc)**2. *Fvm
        elif (v>vc) and (v<=va):
            Fnu=(vm/vc)**(-1/2.) *(va/vm)**(-p/2.) *(v/va)**2.5 *Fvm
        elif (v>va):
            Fnu=(vm/vc)**(-1/2.)*(v/vm)**(-p/2.)*Fvm
    return Fnu

#Fv in slow cooling case        
def Fvslow_sbpl(v,va,vm,vc,Fvm,p):
    if (va<=vm):
        if (v<=va):
            Fnu=(va/vm)**(1/3.)*(v/va)**2. *Fvm
        elif (v>va) and (v<=vm):
            Fnu=(v/vm)**(1/3.)* Fvm
        elif (v>vm) and (v<=vc):
            Fnu=(v/vm)**(-(p-1.)/2.)*Fvm
        elif (v>vc):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    elif (va>vm) and (va<=vc):
        if (v<=vm):
            Fnu=(vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(va/vm)**(-(p-1.)/2.) *(v/va)**2.5 *Fvm
        elif (v>va) and (v<=vc):
            Fnu=(v/vm)**(-(p-1.)/2.)*Fvm
        elif (v>vc):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    elif (va>vc):
        if (v<=vm):
            Fnu=(vc/va)**0.5* (vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(vc/vm)**(-(p-1.)/2.) *(va/vc)**(-p/2.)*(v/va)**2.5 *Fvm
        elif (v>va):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    return Fnu

def Fv_sbpl(v,va,vm,vc,Fvm,p):
    if (vc<=vm):
#        Fnu=Fvfast(v,va,vm,vc,Fvm,p)
        Fnu=Fvfast_sbpl(v,va,vm,vc,Fvm,p)
    else:
#        Fnu=Fvslow(v,va,vm,vc,Fvm,p)
        Fnu=Fvslow_sbpl(v,va,vm,vc,Fvm,p)
    return Fnu



def vFv_sbpl(v,va,vm,vc,Fvm,p):
    return v*Fv_sbpl(v,va,vm,vc,Fvm,p)



#------------------------Simple braken power law for RS---------------------------
#Fv in fast cooling case
def Fvfast_RS(v,va,vm,vc,Fvm,p):
    if (v>vc):
        Fnu=cgs.zero
#        Fnu=np.exp(-100)
    elif (va<=vc):
        if (v<=va):
            Fnu=(va/vc)**(1/3.) *(v/va)**2. *Fvm
        elif (v>va) and (v<=vc):
            Fnu=(v/vc)**(1/3.)* Fvm
        elif (v>vc) and (v<=vm):
            Fnu=(v/vc)**(-1/2.)*Fvm
        elif (v>vm):
            Fnu=(vm/vc)**(-1/2.)*(v/vm)**(-p/2.)*Fvm
    elif (va>vc) and (va<=vm):
        if (v<=vc):
            Fnu=(vc/va)**3. *(v/vc)**2. *Fvm
        elif (v>vc) and (v<=va):
            Fnu=(va/vc)**(-1/2.) *(v/va)**2.5 *Fvm
        elif (v>va) and (v<=vm):
            Fnu=(v/vc)**(-1/2.)*Fvm
        elif (v>vm):
            Fnu=(vm/vc)**(-1/2.)*(v/vm)**(-p/2.)*Fvm    
    elif (va>vm):
        if (v<=vc):
            Fnu=(vc/va)**3. *(vm/va)**((p-1.)/2.) *(v/vc)**2. *Fvm
        elif (v>vc) and (v<=va):
            Fnu=(vm/vc)**(-1/2.) *(va/vm)**(-p/2.) *(v/va)**2.5 *Fvm
        elif (v>va):
            Fnu=(vm/vc)**(-1/2.)*(v/vm)**(-p/2.)*Fvm
    return Fnu

#Fv in slow cooling case        
def Fvslow_RS(v,va,vm,vc,Fvm,p):
    if (v>vc):
        Fnu=cgs.zero
#        Fnu=np.exp(-100)
    elif (va<=vm):
        if (v<=va):
            Fnu=(va/vm)**(1/3.)*(v/va)**2. *Fvm
        elif (v>va) and (v<=vm):
            Fnu=(v/vm)**(1/3.)* Fvm
        elif (v>vm) and (v<=vc):
            Fnu=(v/vm)**(-(p-1.)/2.)*Fvm
#    elif (v>vc):
#                Fnu=np.exp(-700)
#        Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    elif (va>vm) and (va<=vc):
        if (v<=vm):
            Fnu=(vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(va/vm)**(-(p-1.)/2.) *(v/va)**2.5 *Fvm
        elif (v>va) and (v<=vc):
            Fnu=(v/vm)**(-(p-1.)/2.)*Fvm
#    elif (v>vc):
#        Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    elif (va>vc):
#        Fnu=np.exp(-100)
        if (v<=vm):
            Fnu=(vc/va)**0.5* (vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(vc/vm)**(-(p-1.)/2.) *(va/vc)**(-p/2.)*(v/va)**2.5 *Fvm
#    elif (v>va):
#        Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    return Fnu


def Fv_RS(v,va,vm,vc,Fvm,p):
    if (vc<=vm):
        Fnu=cgs.zero
#        Fnu=Fvfast(v,va,vm,vc,Fvm,p)
#        Fnu=Fvfast_RS(v,va,vm,vc,Fvm,p)
    else:
#        Fnu=Fvslow(v,va,vm,vc,Fvm,p)
        Fnu=Fvslow_RS(v,va,vm,vc,Fvm,p)
    return Fnu



def vFv_RS(v,va,vm,vc,Fvm,p):
    return v*Fv_RS(v,va,vm,vc,Fvm,p)




#---------------------Starting update 1 for Fv---------------------------------    
#Current version:
#updates ver.1 for Fv, smooth for nu_a   
#modified for smooth connection between transition at nu_a
#Fv in fast cooling case
def Fvfast(v,va,vm,vc,Fvm,p):
    if (va<=vc):
        tauv = (v/va)**(-5./3.)
        tau = np.exp(-tauv)
        Sv = (va/vc)**(1/3.) *(v/va)**2. *Fvm
        Fnu_thick = Sv*(1.-tau)
        if (tauv <1.e-3):
            Fnu_thick=(v/vc)**(1/3.)* Fvm
        if (v<=va):
                Fnu = Fnu_thick
        elif (v>va) and (v<=vc):
            Fnu = Fnu_thick
        elif (v>vc) and (v<=vm):
            Fnu = Fnu_thick* (v/vc)**(-1/2.-1/3.)
        elif (v>vm):
            Fnu = Fnu_thick* (v/vc)**(-1/2.-1/3.)*(v/vm)**(-p/2.+1/2.)
    elif (va>vc) and (va<=vm):
        tauv = (v/va)**(-3.)
        tau = np.exp(-tauv)
        Sv = (va/vc)**(-1/2.) *(v/va)**2.5 *Fvm
        Fnu_thick = Sv*(1.-tau)
        if (tauv <1.e-3):
            Fnu_thick = (v/vc)**(-1/2.)*Fvm
    if (v<=vc):
        Fnu=Fnu_thick* (v/vc)**(-1/2.)
    elif (v>vc) and (v<=va):
        Fnu = Fnu_thick
    elif (v>va) and (v<=vm):
        Fnu=Fnu_thick
    elif (v>vm):
        Fnu=Fnu_thick* (v/vm)**(-p/2.+1/2.)   
    elif (va>vm):
        tauv = (v/va)**(-(p+5.)/2.)
        tau = np.exp(-tauv)
        Sv = (vm/vc)**(-1/2.) *(va/vm)**(-p/2.) *(v/va)**2.5 *Fvm
        Fnu_thick = Sv*(1.-tau)
        if (tauv <1.e-3):
            Fnu_thick = (vm/vc)**(-1/2.)*(v/vm)**(-p/2.)*Fvm
    if (v<=vc):
        Fnu=Fnu_thick *(v/vc)**(-1/2.)
    elif (v>vc) and (v<=va):
        Fnu=Fnu_thick
    elif (v>va):
        Fnu=Fnu_thick
    return Fnu

#Fv in slow cooling case        
def Fvslow(v,va,vm,vc,Fvm,p):
    if (va<=vm):
        tauv = (v/va)**(-5./3.)
        tau = np.exp(-tauv)
        Sv = (va/vm)**(1/3.)*(v/va)**2. *Fvm
        Fnu_thick = Sv*(1.-tau)
        if (tauv <1.e-3):
            Fnu_thick = (v/vm)**(1/3.)* Fvm
        if (v<=va):
            Fnu= Fnu_thick
        elif (v>va) and (v<=vm):
            Fnu= Fnu_thick
        elif (v>vm) and (v<=vc):
            Fnu= Fnu_thick* (v/vm)**(-(p-1.)/2.-1/3.)
        elif (v>vc):
            Fnu=Fnu_thick* (v/vm)**(-(p-1.)/2.-1/3.)*(v/vc)**(-1/2.)
    elif (va>vm) and (va<=vc):
        tauv = (v/va)**(-(p+4.)/2.)
        tau = np.exp(-tauv)
        Sv = (va/vm)**(-(p-1.)/2.) *(v/va)**2.5 *Fvm
        Fnu_thick = Sv*(1.-tau)
        if (tauv <1.e-3):
            Fnu_thick = (v/vm)**(-(p-1.)/2.)*Fvm
        if (v<=vm):
            Fnu= Fnu_thick *(v/vm)**(-1/2.)
        elif (v>vm) and (v<=va):
            Fnu= Fnu_thick
        elif (v>va) and (v<=vc):
            Fnu= Fnu_thick
        elif (v>vc):
            Fnu= Fnu_thick* (v/vc)**(-1/2.)
    elif (va>vc):
        tauv = (v/va)**(-(p+5.)/2.)
        tau = np.exp(-tauv)
        Sv = (vc/vm)**(-(p-1.)/2.) *(va/vc)**(-p/2.)*(v/va)**2.5 *Fvm
        Fnu_thick = Sv*(1.-tau)
        if (tauv <1.e-3):
            Fnu_thick = (vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
        if (v<=vm):
            Fnu= Fnu_thick* (v/vm)**(-1/2.)
        elif (v>vm) and (v<=va):
            Fnu= Fnu_thick
        elif (v>va):
            Fnu= Fnu_thick
    return Fnu

#--------------------------end of update 1 for Fv------------------------------




#--------------------------Start update 2 for Fv------------------------------
#smooth for vm, vc:
#status: ongoing

def Fvslow2(v,va,vm,vc,Fvm,p):
    if (va<=vm):
        if (v<=va):
            Fnu=(va/vm)**(1/3.)*(v/va)**2. *Fvm
        elif (v>va):
            p1=p
            p2=p + 1.
            x1 = v/vm
            x2 = v/vc
            if (x2 < 10.): 
                if (x1 >20.):    x1 = 20.
    #                    s1= math.gamma(p1/4.+19/12.)*math.gamma(p1/4.-1/12.)
                f1 = lambda x: Fx(x) * x**((p1-3.)/2.)
                s1 = quad(f1,x2,x1)[0]
            else:
                s1 = 0.
                if (x2 > 20.): x2=20.
    #                s2= math.gamma(p2/4.+19/12.)*math.gamma(p2/4.-1/12.)
            f2 = lambda x: Fx(x) * x**((p2-3.)/2.)
            s2 = quad(f2,x2/100.,x2)[0]

                
            Fnu1 = (v/vm)**(-(p1-1.)/2.) * s1    
            Fnu2 =(vc/vm)**(-(p1-1.)/2.)*(v/vc)**(-(p2-1.)/2.)*s2

            Fnu = (Fnu1 + Fnu2) * Fvm
    elif (va>vm) and (va<=vc):
        if (v<=vm):
            Fnu=(vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(va/vm)**(-(p-1.)/2.) *(v/va)**2.5 *Fvm
        elif (v>va) and (v<=vc):
            Fnu=(v/vm)**(-(p-1.)/2.)*Fvm
        elif (v>vc):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    elif (va>vc):
        if (v<=vm):
            Fnu=(vc/va)**0.5* (vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(vc/vm)**(-(p-1.)/2.) *(va/vc)**(-p/2.)*(v/va)**2.5 *Fvm
        elif (v>va):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    return Fnu

#the third version
#modified for smooth connection between transition        
def Fvslow3(v,va,vm,vc,Fvm,p):
    if (va<=vm):
        if (v<=va):
            Fnu=(va/vm)**(1/3.)*(v/va)**2. *Fvm
        elif (v>va) and (v<=vm):
            p1=p
            p2=p + 1.
            x1 = v/vm
            x2 = v/vc
            f1 = lambda x: Fx(x) * x**((p1-3.)/2.)
            s = quad(f1,x2,x1)[0]
            Fnu1 = (v/vm)**(-(p1-1.)/2.) * s

            f2 = lambda x: Fx(x) * x**((p2-3.)/2.)
            s = quad(f2,x2/10.,x2)[0]
            Fnu2 =(vc/vm)**(-(p1-1.)/2.)*(v/vc)**(-(p2-1.)/2.)*s

            Fnu = (Fnu1 + Fnu2) * Fvm
        elif (v>vm) and (v<=vc):
            p1=p
            p2=p + 1.
            x1 = 1.
            x2 = v/vc
            f1 = lambda x: Fx(x) * x**((p1-3.)/2.)
            s = quad(f1,x2,x1)[0]
            Fnu1 = (v/vm)**(-(p1-1.)/2.) * s

            f2 = lambda x: Fx(x) * x**((p2-3.)/2.)
            s = quad(f2,x2/10.,x2)[0]
            Fnu2 =(vc/vm)**(-(p1-1.)/2.)*(v/vc)**(-(p2-1.)/2.)*s

            Fnu = (Fnu1 + Fnu2) * Fvm
        elif (v>vc):
            p1=p
            p2=p + 1.
            f2 = lambda x: Fx(x) * x**((p2-3.)/2.)
            s = quad(f2,0.01,1.)[0]
            Fnu2 =(vc/vm)**(-(p1-1.)/2.)*(v/vc)**(-(p2-1.)/2.)*s

            Fnu = (Fnu2) * Fvm
    elif (va>vm) and (va<=vc):
        if (v<=vm):
            Fnu=(vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(va/vm)**(-(p-1.)/2.) *(v/va)**2.5 *Fvm
        elif (v>va) and (v<=vc):
            Fnu=(v/vm)**(-(p-1.)/2.)*Fvm
        elif (v>vc):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    elif (va>vc):
        if (v<=vm):
            Fnu=(vc/va)**0.5* (vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(vc/vm)**(-(p-1.)/2.) *(va/vc)**(-p/2.)*(v/va)**2.5 *Fvm
        elif (v>va):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    return Fnu

def Get_N(a,b,width):
    N=int((b-a)/width + 1)
    if N%2 == 0:
        N=N+1
    return N


def Fvslow4(v,va,vm,vc,Fvm,p):
    if (va<=vm):
        if (v<=va):
            Fnu=(va/vm)**(1/3.)*(v/va)**2. *Fvm
        elif (v>va):
            p1=p
            p2=p + 1.
            x1 = v/vm
            x2 = v/vc
            f1 = lambda x: Fx(x) * x**((p1-3.)/2.)
            datas = []
            r=x2
            width = (x1-x2)/10.
            inum=Get_N(x2,x1,width)
            for i in range(0,inum):
                datas.append(f1(r))
                r = r+width
            sum = datas[0]+datas[inum-1]
            for i in range(2,inum):
                if i%2== 0:
                    sum = sum +4*datas[i-1]
                else:
                    sum = sum +2*datas[i-1]
            s=sum*width/3.0
            Fnu1 = (v/vm)**(-(p1-1.)/2.) * s

            f2 = lambda x: Fx(x) * x**((p2-3.)/2.)
            datas = []
            r=0.
            width = (x2-0.)/10.
            inum=Get_N(0.,x2,width)
            for i in range(0,inum):
                datas.append(f2(r))
                r = r+width
            sum = datas[0]+datas[inum-1]
            for i in range(2,inum):
                if i%2== 0:
                    sum = sum +4*datas[i-1]
                else:
                    sum = sum +2*datas[i-1]
            s=sum*width/3.0
            Fnu2 =(vc/vm)**(-(p1-1.)/2.)*(v/vc)**(-(p2-1.)/2.)*s

            Fnu = (Fnu1 + Fnu2) * Fvm
    elif (va>vm) and (va<=vc):
        if (v<=vm):
            Fnu=(vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(va/vm)**(-(p-1.)/2.) *(v/va)**2.5 *Fvm
        elif (v>va) and (v<=vc):
            Fnu=(v/vm)**(-(p-1.)/2.)*Fvm
        elif (v>vc):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    elif (va>vc):
        if (v<=vm):
            Fnu=(vc/va)**0.5* (vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(vc/vm)**(-(p-1.)/2.) *(va/vc)**(-p/2.)*(v/va)**2.5 *Fvm
        elif (v>va):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    return Fnu

def Fvslow5(v,va,vm,vc,Fvm,p):
    if (va<=vm):
        if (v<=va):
            Fnu=(va/vm)**(1/3.)*(v/va)**2. *Fvm
        elif (v>va):
            p1=p
            p2=p + 1.
            x1 = v/vm
            x2 = v/vc
            if (x2 < 5.): 
                if (x1 >10.): 
                    s1= math.gamma(p1/4.+19/12.)*math.gamma(p1/4.-1/12.)
                else:
                    f1 = lambda x: Fx(x) * x**((p1-3.)/2.)
                    s1 = quad(f1,x2,x1)[0]
                if (x2 < 0.1): 
                    s2 = 0.
                else:
                    f2 = lambda x: Fx(x) * x**((p2-3.)/2.)
                    s2 = quad(f2,x2/100.,x2)[0]
            else:
                s1 = 0.
                s2= math.gamma(p2/4.+19/12.)*math.gamma(p2/4.-1/12.)
                
            Fnu1 = (v/vm)**(-(p1-1.)/2.) * s1    
            Fnu2 =(vc/vm)**(-(p1-1.)/2.)*(v/vc)**(-(p2-1.)/2.)*s2

            Fnu = (Fnu1 + Fnu2) * Fvm
    elif (va>vm) and (va<=vc):
        if (v<=vm):
            Fnu=(vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(va/vm)**(-(p-1.)/2.) *(v/va)**2.5 *Fvm
        elif (v>va) and (v<=vc):
            Fnu=(v/vm)**(-(p-1.)/2.)*Fvm
        elif (v>vc):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    elif (va>vc):
        if (v<=vm):
            Fnu=(vc/va)**0.5* (vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(vc/vm)**(-(p-1.)/2.) *(va/vc)**(-p/2.)*(v/va)**2.5 *Fvm
        elif (v>va):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    return Fnu


def Fvslow6(v,va,vm,vc,Fvm,p):
    if (va<=vm):
        if (v<=va):
            Fnu=(va/vm)**(1/3.)*(v/va)**2. *Fvm
        elif (v>va):
            p1=p
            p2=p + 1.
            x1 = v/vm
            x2 = v/vc
            if (x1 <= 0.1):
                Fnu1 = (v/vm)**(1/3.)
                s2 =0.
            elif (x1 < 5.):
                f1 = lambda x: Fx(x) * x**((p1-3.)/2.)
                s1 = quad(f1,x2,x1)[0]
                Fnu1 = (v/vm)**(-(p1-1.)/2.)*s1
            else:
                Fnu1 = (v/vm)**(-(p1-1.)/2.)
                
            if (x2 > 5.):
                s2 =1.
            elif (x1 > 0.1):
                f2 = lambda x: Fx(x) * x**((p2-3.)/2.)
                s2 = quad(f2,x2/100.,x2)[0]

            Fnu2 =(vc/vm)**(-(p1-1.)/2.)*(v/vc)**(-(p2-1.)/2.)*s2

            Fnu = (Fnu1 + Fnu2) * Fvm
    elif (va>vm) and (va<=vc):
        if (v<=vm):
            Fnu=(vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(va/vm)**(-(p-1.)/2.) *(v/va)**2.5 *Fvm
        elif (v>va) and (v<=vc):
            Fnu=(v/vm)**(-(p-1.)/2.)*Fvm
        elif (v>vc):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    elif (va>vc):
        if (v<=vm):
            Fnu=(vc/va)**0.5* (vm/va)**((p+4.)/2.)*(v/vm)**2. * Fvm
        elif (v>vm) and (v<=va):
            Fnu=(vc/vm)**(-(p-1.)/2.) *(va/vc)**(-p/2.)*(v/va)**2.5 *Fvm
        elif (v>va):
            Fnu=(vc/vm)**(-(p-1.)/2.)*(v/vc)**(-p/2.)*Fvm
    return Fnu

def Fv1(v,va,vm,vc,Fvm,p):
    if (vc<=vm):
        Fnu=Fvfast1(v,va,vm,vc,Fvm,p)
    else:
        Fnu=Fvslow1(v,va,vm,vc,Fvm,p)
    return Fnu

def Fv(v,va,vm,vc,Fvm,p):
    if (vc<=vm):
        Fnu=Fvfast(v,va,vm,vc,Fvm,p)
#        Fnu=Fvfast_sbpl(v,va,vm,vc,Fvm,p)
    else:
        Fnu=Fvslow(v,va,vm,vc,Fvm,p)
#        Fnu=Fvslow_sbpl(v,va,vm,vc,Fvm,p)
    return Fnu



def vFv(v,va,vm,vc,Fvm,p):
    return v*Fv(v,va,vm,vc,Fvm,p)

#--------------------------------------------------------------------------
#---------------------------------SSC--------------------------------------
#--------------------------------------------------------------------------
#-----------------------Trail version of SSC
def FvIC_trail(v, gm_m, gm_c,vm, vc, Fvm_IC, p):
    if (vc<=vm):
#------fast cooling-------
        vm_IC=gm_c**2 *vm
        vc_IC=gm_m**2 *vc
        if (v<=vc_IC):
            Fnu=(v/vc_IC)**(1/3.)* Fvm_IC
        elif (v>vc_IC) and (v<=vm_IC):
            Fnu=(v/vc_IC)**(-1/2.)* Fvm_IC
        elif (v>vm_IC):
            Fnu=(vm_IC/vc_IC)**(-1/2.)*(v/vm_IC)**(-p/2.)*Fvm_IC        
    else:
#------slow cooling------
        vm_IC=gm_m**2 *vm
        vc_IC=gm_c**2 *vc
        if (v<=vm_IC):
            Fnu=(v/vm_IC)**(1/3.)* Fvm_IC
        elif (v>vm_IC) and (v<=vc_IC):
            Fnu=(v/vm_IC)**(-(p-1.)/2.)* Fvm_IC
        elif (v>vc_IC):
            Fnu=(vc_IC/vm_IC)**(-(p-1.)/2.)*(v/vc_IC)**(-p/2.)* Fvm_IC
    return Fnu

#-----------------------Simple broken power law (sbpl) version of SSC
#Fv in fast cooling case
def FvICfast_sbpl(v, gm_a,gm_m, gm_c,va,vm, vc,Fvm_IC, p):
    x0=0.4714
    vca=4.*gm_c**2 *va *x0
    vcm=4.*gm_c**2 *vm *x0
    vcc=4.*gm_c**2 *vc *x0
    vmm=4.*gm_m**2 *vm *x0
    vam=4.*gm_a**2 *vm *x0
    vaa=4.*gm_a**2 *va *x0
#    vmmax=4.*gm_m**2 *vm *x0
    if (va<=vc):
        if (v<=vca):
            Fnu=(vca/vcc)**(1/3.) *(v/vca)
        elif (v>vca) and (v<=vcc):
            Fnu=(v/vcc)**(1/3.)
        elif (v>vcc) and (v<=vmm):
            Fnu=(v/vcc)**(-1/2.)
        elif (v>vmm):
            Fnu=(vmm/vcc)**(-1/2.)*(v/vmm)**(-p/2.)
    elif (va>vc) and (va<=vm):
        if (v<=vaa):
            Fnu=(v/vaa)
        elif (v>vaa) and (v<=vmm):
            Fnu=(v/vaa)**(-1/2.)
        elif (v>vmm):
            Fnu=(vmm/vaa)**(-1/2.)*(v/vmm)**(-p/2.)    
    elif (va>vm):
        if (v<=vaa):
            Fnu=(v/vaa)
        else:
            Fnu=(v/vaa)**(-p/2.)
    return Fnu* Fvm_IC


#Fv in slow cooling case        
def FvICslow_sbpl(v, gm_a,gm_m, gm_c, va, vm, vc, Fvm_IC, p):
    x0=0.4714
    vca=4.*gm_c**2 *va *x0
    vma=4.*gm_m**2 *va *x0
    vmm=4.*gm_m**2 *vm *x0
    vmc=4.*gm_m**2 *vc *x0
    vaa=4.*gm_a**2 *va *x0
    vcc=4.*gm_c**2 *vc *x0
    if (va<=vm):
        if (v<=vma):
            Fnu=(vma/vmm)**(1/3.)*(v/vma) 
        elif (v>vma) and (v<=vmm):
            Fnu=(v/vmm)**(1/3.)
        elif (v>vmm) and (v<=vcc):
            Fnu=(v/vmm)**(-(p-1.)/2.)
        elif (v>vcc):
            Fnu=(vcc/vmm)**(-(p-1.)/2.)*(v/vcc)**(-p/2.)
    elif (va>vm) and (va<=vc):
        if (v<=vma):
            Fnu=(vmm/vma)**((p+1.)/2.)*(v/vmm)
        elif (v>vma) and (v<=vcc):
            Fnu=(v/vmm)**(-(p-1.)/2.)
        elif (v>vcc):
            Fnu=(vcc/vmm)**(1/2.)*(v/vmm)**(-p/2.)
    elif (va>vc):
        if (v<=vaa):
            Fnu=(v/vaa)
        else:
            Fnu=(v/vaa)**(-p/2.)
    return Fnu* Fvm_IC


#------------------------He & Lei's version of SSC------------
#Fv in fast cooling case
def FvICfast(v, gm_a,gm_m, gm_c,va,vm, vc, Fvm_IC, p):
    x0=0.4714
    vca=4.*gm_c**2 *va *x0
    vcm=4.*gm_c**2 *vm *x0
    vcc=4.*gm_c**2 *vc *x0
    vmm=4.*gm_m**2 *vm *x0
    vam=4.*gm_a**2 *vm *x0
    vaa=4.*gm_a**2 *va *x0
    if (va<=vc):
        if (v<=vca):
            Fnu=(va/vc)**(1/3.) *(v/vca) *Fvm_IC*5/6.
        elif (v>vca) and (v<=vcc):
            Fnu=(v/vcc)**(1/3.)* Fvm_IC*9/10.
        elif (v>vcc) and (v<=vcm):
            Fnu=(v/vcc)**(-1/2.)*Fvm_IC*(28/15.+np.log(v/vcc))/3.
        elif (v>vcm) and (v<=vmm):
            Fnu=(v/vcc)**(-1/2.)*Fvm_IC*(2*(p+5.)/(p+2.)/(p-1.)-2*(p-1.)/3./(p+2.)+np.log(vmm/v))/3.
        elif (v>vmm):
            Fnu=(vc/vm)*(v/vmm)**(-p/2.)*Fvm_IC*(2*(p+5.)/(p-1.)/3.-2*(p-1.)/(p+2.)/3.+np.log(v/vmm))/(p+2.)
    elif (va>vc) and (va<=vm):
        R=gm_c/gm_a/3.
        if (v<=vaa):
            Fnu=(0.5*R+1.)*(R+4.)*(v/vaa) *Fvm_IC
        elif (v>vaa) and (v<=vam):
            Fnu=R*(v/vaa)**(-1/2.)*Fvm_IC*(R/6.+0.9+R*np.log(v/vaa)/4.)
        elif (v>vam) and (v<=vmm):
            Fnu=R**2.*(v/vaa)**(-1/2.)*Fvm_IC*(3/(p-1.)-0.5+3.*np.log(vmm/v)/4.)
        elif (v>vmm):
            Fnu=R**2.*(va/vm)*(v/vmm)**(-p/2.)*Fvm_IC*(4/(p+3.)*(gm_a/gm_m)**(p-1.)*gm_a/gm_c+3*(p+1.)/(p-1./(p-2.+0.5*np.log(v/vmm))))*9./2./(p+2.)    
    elif (va>vm):
        R=(gm_m/gm_a)**(p-1.)*gm_c/gm_a/3.
        if (v<=vaa):
            Fnu=(3*R/2./(p+2.)+1.)*(3*R/(p+2.)+4.)*(v/vaa) *Fvm_IC
        elif (v>vaa):
            Fnu=(6*R/(p+3.)+R*(9*R/2/(p+2.)+1.)+9*R/4.*np.log(v/vaa))*(v/vaa)**(-p/2.)*Fvm_IC/(p+2.)
    return Fnu

#Fv in slow cooling case        
def FvICslow(v, gm_a,gm_m, gm_c, va, vm, vc, Fvm_IC, p):
    x0=0.4714
    vca=4.*gm_c**2 *va *x0
    vma=4.*gm_m**2 *va *x0
    vmm=4.*gm_m**2 *vm *x0
    vmc=4.*gm_m**2 *vc *x0
    vaa=4.*gm_a**2 *va *x0
    vcc=4.*gm_c**2 *vc *x0
    if (va<=vm):
        if (v<=vma):
            Fnu=(va/vm)**(1/3.)*(v/vma) *Fvm_IC*5*(p-1.)/(2*(p+1.))
        elif (v>vma) and (v<=vmm):
            Fnu=(v/vmm)**(1/3.)* Fvm_IC*3*(p-1.)/(2*(p-1/3.))
        elif (v>vmm) and (v<=vmc):
            Fnu=(v/vmm)**(-(p-1.)/2.)*Fvm_IC*(4*(p+1/3.)/(p+1.)/(p-1/3.)+np.log(v/vmm))*(p-1.)/(p+1.)
        elif (v>vmc) and (v<=vcc):
            Fnu=(v/vmm)**(-(p-1.)/2.)*Fvm_IC*(2*(2*p+3.)/(p+2.)-2/(p+1.)/(p+2.)+np.log(vcc/v))*(p-1.)/(p+1.)
        elif (v>vcc):
            Fnu=(v/vmm)**(-p/2.)*(vc/vm)*Fvm_IC*(2*(2*p+3.)/(p+2.)-2/(p+2.)**2.+np.log(v/vcc)*(p+1.)/(p+2.))*(p-1.)/(p+1.)
#        Fnu=(vcc/vmm)**(-(p-1.)/2.)*(v/vcc)**(-p/2.)*Fvm_IC
    elif (va>vm) and (va<=vc):
        if (v<=vma):
            Fnu=(vm/va)**((p+1.)/2.)*(v/vmm) * Fvm_IC*2*(p+4.)*(p-1.)/3./(p+1.)**2.
        elif (v>vma) and (v<=vmc):
            Fnu=(v/vmm)**(-(p-1.)/2.)*Fvm_IC*(2*(2*p+5.)/(p+1.)/(p+4.)+ np.log(v/vma))*(p-1.)/(p+1.)
        elif (v>vmc) and (v<=vca):
            Fnu=(v/vmm)**(-(p-1.)/2.)*Fvm_IC*(2.+2./(p+4.)+np.log(vc/va))*(p-1.)/(p+1.)
        elif (v>vca) and (v<=vcc):
            Fnu=(v/vmm)**(-(p-1.)/2.)*Fvm_IC*(2*(2*p+1.)/(p+1.)+np.log(vcc/v))*(p-1.)/(p+1.)
        elif (v>vcc):
            Fnu=(vc/vm)*(v/vmm)**(-p/2.)*Fvm_IC*(2*(2*p+5.)/(p+2.)+np.log(v/vcc))*(p-1.)/(p+1.)
    elif (va>vc):
        R=(p-1.)*(gm_m/gm_a)**(p-1.)*gm_c/gm_a/3.
        if (v<=vaa):
            Fnu=(3*R/2./(p+2.)+1.)*(3*R/(p+2.)+4.)*(v/vaa) *Fvm_IC
        elif (v>vaa):
            Fnu=(6*R/(p+3.)+R*(9*R/2/(p+2.)+1.)+9*R/4.*np.log(v/vaa))*(v/vaa)**(-p/2.)*Fvm_IC/(p+2.)
    return Fnu

#no KN, no vMax effects
def FvIC0(v, gm_a,gm_m, gm_c, va, vm, vc, Fvm_IC, p):
    if (vc<=vm):
#        Fnu=FvICfast(v, gm_a,gm_m, gm_c, va, vm, vc, Fvm_IC, p)
        Fnu=FvICfast_sbpl(v, gm_a,gm_m, gm_c, va, vm, vc, Fvm_IC, p)
        
    else:
#        Fnu=FvICslow(v, gm_a,gm_m, gm_c, va, vm, vc, Fvm_IC, p)
        Fnu=FvICslow_sbpl(v, gm_a,gm_m, gm_c, va, vm, vc, Fvm_IC, p)
    return Fnu


def FvIC(v, gm_a,gm_m, gm_c, va, vm, vc, vMax, Fvm_IC, GM, p):
    if (vc<=vm):
#        Fnu=FvICfast(v, gm_a,gm_m, gm_c, va, vm, vc, Fvm_IC, p)
        Fnu=FvICfast(v, gm_a,gm_m, gm_c, va, vm, vc, Fvm_IC, p)
        x=cgs.h*v/(GM*cgs.me*cgs.c**2. *gm_m)
    else:
#        Fnu=FvICslow(v, gm_a,gm_m, gm_c, va, vm, vc, Fvm_IC, p)
        Fnu=FvICslow(v, gm_a,gm_m, gm_c, va, vm, vc, Fvm_IC, p)
        x=cgs.h*v/(GM*cgs.me*cgs.c**2. *gm_c)
    sigmaR=astro.sigmaKN(x)/cgs.sigmaT
    return Fnu*sigmaR
