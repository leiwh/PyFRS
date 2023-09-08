#==================================Physical Constants useful in Astronomy, in cgs units.=================================
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

#!/usr/bin/python
# Filename: cgs.py
from astropy import constants as const
import numpy as np

def about():
	print("List of the constants in cgs")
version = '2.0'

def help():
    print("list")

zero=np.exp(-700)
day=24.0*3600.0
year=day*365.0
deg=180.0/3.14159
pi=3.141592653589793

#Speed of light
#c=2.99792458*1.e10
c=const.c.cgs.value
#Plank's constant
#h= 6.6260755e-27
h=const.h.cgs.value

#Reduced Planck constant
#hbar=h/2pi
#hbar=1.05457266e-27
hbar=const.hbar.cgs.value

#Boltzmann's constant
#k=1.380658e-16
k_B=const.k_B.cgs.value
k=k_B

#Electron charge
#eq=4.8e-10
e=const.e.gauss.value
eq=e
#qe=4.8e-10
qe=eq

#---mass
#electron mass
#me=9.11e-28
m_e=const.m_e.cgs.value
me=m_e
#proton mass
#mp=1.6726231e-24
mp=const.m_p.cgs.value
#neutron mass
#mn=1.674929e-24
mn=const.m_n.cgs.value
mH=1.673534e-24

#Gravitational constant
#G=6.67259e-8
G=const.G.cgs.value

#Avogadro number
#NA=6.022e23
N_A=const.N_A.cgs.value
NA=N_A

#Atomic mass
u=1.66e-24
#u=const.u.cgs.value

#Fine-structure constant
aa=7.297e-3
alpha= aa
#alpha=const.alpha.cgs.value
aa=alpha

qm=5.273e17
R0=1.097e5
a0=5.29e-9
lamda=2.43e-10
re=2.82e-13

#Gas constant
#R=8.31e7
R=const.R.cgs.value

#Stefan-Boltzmann constant
#sigma=5.67e-5
sigma_sb=const.sigma_sb.cgs.value
sigma=sigma_sb

#Thomson cross section
sigmaT=6.65e-25
#Radiation constant: a=4sigma/c
a=7.566e-15


#---distance
#Astronomical Unit 
#AU=1.496e13
au=const.au.cgs.value
AU=au
#A=1.496e13
A=au

#Kiloparsec
#pc=3.086e18
pc=const.pc.cgs.value
#kpc=3.086e21
kpc=const.kpc.cgs.value
Mpc=3.085677581467192e24

#pc2A=206264.806
#pc2AU=206264.806
pc2AU=const.pc.to('au').value
pc2A=pc2AU

pc2ly=3.26
ly=9.461e17
ly2A=6.324e4
ly2pc=0.3066
mile=160934.4

#--energy, frequence, flux
eV2cm=12398.52e-8
eV2bs=8065.479
eV2Hz=2.418e14
keV2Hz=2.418e17
MeV2Hz=2.418e20
eV2erg=1.6e-12
keV2erg=1.6e-9
MeV2erg=1.6e-6
eV2K=11604.5
keV2K=eV2K*1e3
MeV2K=eV2K*1e6
K2eV=1./eV2K
K2keV=1./keV2K
K2MeV=1./MeV2K
erg=1.e-7
esu=0.33e-9
dyn=1.e-5
mecc2MeV=0.511
Jy=1.e-23
mJy=1.e-26
uJy=1.e-29
Crab=2.4e-8
mCrab=2.4e-11


#---Sun, Earth, Jupiter
#Solar mass
#Msun=1.9891e33
M_sun=const.M_sun.cgs.value
Msun=M_sun

#Solar radius
#Rsun=6.9599e10
R_sun=const.R_sun.cgs.value
Rsun=R_sun

#Solar distance to earth
Dsun2earth=1.5e13

#Solar luminosity
#Lsun=3.826e33
L_sun=const.L_sun.cgs.value
Lsun=L_sun

#distance from Sun to Sgr A*
Dsun=8.33*kpc

#Earth mass
#Mearth=6.0e27
M_earth=const.M_earth.cgs.value
Mearth=M_earth

#Earth radius
#Rearth=6.4e8
R_earth=const.R_earth.cgs.value
Rearth=R_earth

#Jupiter mass
#MJupiter=1.9e30
M_jup=const.M_jup.cgs.value
MJupiter=M_jup

#Jupiter radius
#RJupiter=7.1e9
R_jup=const.R_jup.cgs.value
RJupiter=R_jup

#Moon mass
Mmoon=7.4e25
#Moon radius
Rmoon=1.7e8
#Moon distance to earth
Dmoon=3.8e10

#-------other
#inch to cm
inch=2.539999918


#---Other unit in astrophysics
Msun2erg=Msun* c**2.
rg=G* Msun/c**2.
UnitOmega=c**3/G/Msun


# End of units.py
