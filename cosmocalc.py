#!/usr/bin/env python
"""
Calculate useful values for a given cosmology.  This module uses code adapted
from `CC.py`_ (`James Schombert`_) which is a Python version of the
`Cosmology Calculator`_ (`Ned Wright`_).

The following values are calculated:

    ====  ===================================  ===========
    Name  Value                                Units
    ====  ===================================  ===========
    z     Input redshift 
    H0    Hubble constant
    WR    Omega(radiation)                     
    WK    Omega curvaturve = 1-Omega(total)     
    WM    Omega matter
    WV    Omega vacuum
    DTT   Time from z to now                   Gyr
    age   Age of Universe                      Gyr
    zage  Age of Universe at redshift z        Gyr
    DCMR  Comoving radial distance             Gyr Mpc cm
    VCM   Comoving volume within redshift      Gpc3 
    DA    Angular size distance                Gyr Mpc cm
    DL    Luminosity distance                  Gyr Mpc cm
    PS    Plate scale - distance per arcsec    kpc cm
    ====  ===================================  ===========

.. _`James Schombert`: http://abyss.uoregon.edu/~js/
.. _`CC.py`: http://www.astro.ucla.edu/~wright/CC.python
.. _`Ned Wright`: http://www.astro.ucla.edu/~wright/intro.html
.. _`Cosmology Calculator`: http://www.astro.ucla.edu/~wright/CosmoCalc.html

:Copyright: Smithsonian Astrophysical Observatory (2009)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

import math

# Define a few constants
cm_per_pc = 3.0856775813057289536e+18 
c = 299792.458                         # velocity of light in km/sec
km_per_ly = 3600*24*365.25*c           # km per light-year
Tyr = 977.8                            # coefficent for converting 1/H into Gyr
arcsec_per_rad = 206264.806
_outvals_str = ('z H0 WM WV WK WR',
                'DA DA_Gyr DA_Mpc DA_cm',
                'DL DL_Gyr DL_Mpc DL_cm',
                'DCMR DCMR_Gyr DCMR_Mpc DCMR_cm',
                'PS_kpc PS_cm',
                'DTT DTT_Gyr',
                'VCM VCM_Gpc3',
                'age age_Gyr',
                'zage zage_Gyr',)
_outvals = (' '.join(_outvals_str)).split()

def cosmocalc(z, H0=71, WM=0.27, WV=None):
    """
    Calculate useful values for the supplied cosmology.

    This routine returns a dictionary of values in the form ``<name>: <value>``,
    where the values are supplied in "natural" units for cosmology, e.g. 1/H0.
    In addition various useful unit conversions are done and stored in the
    dictionary as ``<name>_<unit>: <value>``.  E.g. angular size distance::

      'DA': 0.38250549415474988,
      'DA_Gyr': 5.2678010166833023,
      'DA_Mpc': 1615.1022857909447,
      'DA_cm': 4.9836849147807571e+27

    Example::

     >>> from cosmocalc import cosmocalc
     >>> from pprint import pprint
     >>> pprint(cosmocalc(3, H0=75, WM=.25))
     {'DA': 0.39103776375786625,
      'DA_Gyr': 5.0980896720325548,
      'DA_Mpc': 1563.0689649039205,
      'DA_cm': 4.8231268630387788e+27,
      'DCMR': 1.564151055031465,
      'DCMR_Gyr': 20.392358688130219,
      'DCMR_Mpc': 6252.2758596156818,
      'DCMR_cm': 1.9292507452155115e+28,
      'DL': 6.25660422012586,
      'DL_Gyr': 81.569434752520877,
      'DL_Mpc': 25009.103438462727,
      'DL_cm': 7.717002980862046e+28,
      'DTT': 0.84826379084317027,
      'DTT_Gyr': 11.059097795819358,
      'H0': 75,
      'PS_cm': 2.3383178917293232e+22,
      'PS_kpc': 7.5779721961095019,
      'VCM': 1.2756009121294902,
      'VCM_Gpc3': 1023.7714254161302,
      'WK': 0.0,
      'WM': 0.25,
      'WR': 7.4044444444444448e-05,
      'WV': 0.74992595555555552,
      'age': 1.0133755371756261,
      'age_Gyr': 13.211714670004362,
      'z': 3,
      'zage': 0.16511174633245579,
      'zage_Gyr': 2.1526168741850036}

    :param z: redshift
    :param H0: Hubble constant (default = 71)
    :param WM: Omega matter (default = 0.27)
    :param WV: Omega vacuum (default = 1.0 - WM - 0.4165/(H0*H0))

    :rtype: dictionary of cosmology values (name_unit = value)
    """

    if z > 100:
        z = z / 299792.458    # Values over 100 are in km/s

    if WV is None:
        WV = 1.0 - WM - 0.4165/(H0*H0)  # Omega(vacuum) or lambda
    age = 0.0      # age of Universe in units of 1/H0
    
    h = H0 / 100.
    WR = 4.165E-5 / (h*h)   # includes 3 massless neutrino species, T0 = 2.72528
    WK = 1 - WM - WR - WV
    az = 1.0 / (1 + 1.0*z)
    n=1000         # number of points in integrals
    for i in range(n):
        a = az * (i + 0.5) / n
        adot = math.sqrt(WK + (WM/a) + (WR/(a*a)) + (WV*a*a))
        age = age + 1./adot
    
    zage = az * age / n
    DTT = 0.0
    DCMR = 0.0
    
    # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    for i in range(n):
        a = az + (1-az) * (i+0.5) / n
        adot = math.sqrt(WK + (WM/a) + (WR/(a*a)) + (WV*a*a))
        DTT = DTT + 1./adot
        DCMR = DCMR + 1./(a*adot)
    
    DTT = (1.-az) * DTT / n
    DCMR = (1.-az) * DCMR / n
    age = DTT + zage
    
    # tangential comoving distance
    ratio = 1.0
    x = math.sqrt(abs(WK)) * DCMR
    if x > 0.1:
        if WK > 0:
            ratio =  0.5 * (math.exp(x) - math.exp(-x)) / x 
        else:
            ratio = math.math.sin(x) / x
    else:
        y = x * x
        if WK < 0:
            y = -y
        ratio = 1. + y/6. + y*y/120.
    DCMT = ratio * DCMR

    # comoving volume computation
    ratio = 1.00
    x = math.sqrt(abs(WK)) * DCMR
    if x > 0.1:
        if WK > 0:
            ratio = (0.125 * (math.exp(2.*x) - math.exp(-2.*x)) - x/2.) / (x**3 / 3.)
        else:
            ratio = (x/2. - math.sin(2.*x)/4.)/(x**3 / 3.)
    else:
        y = x * x
        if WK < 0: y = -y
        ratio = 1. + y/5. + (2./105.)*y*y
    VCM = ratio * DCMR**3 / 3.
    VCM_Gpc3 = 4. * math.pi * (0.001*c/H0)**3 * VCM

    DA = az * DCMT
    DL = DA / (az*az)

    # Now convert to some more useful units
    Gyr = lambda x: Tyr / H0 * x
    Mpc = lambda x: c / H0 * x
    cm = lambda x: Mpc(x) * 1e6 * cm_per_pc
    
    DA_Gyr = Gyr(DA)
    DA_Mpc = Mpc(DA)
    DA_cm = cm(DA)

    DL_Gyr = Gyr(DL)
    DL_Mpc = Mpc(DL)
    DL_cm = cm(DL)

    DCMR_Gyr = Gyr(DCMR)
    DCMR_Mpc = Mpc(DCMR)
    DCMR_cm = cm(DCMR)

    DTT_Gyr = Gyr(DTT)
    age_Gyr = Gyr(age)
    zage_Gyr = Gyr(zage)
    
    PS_kpc = Mpc(DA) * 1000 / arcsec_per_rad
    PS_cm = PS_kpc * cm_per_pc * 1000

    localvals = locals()
    return dict((x, localvals[x]) for x in _outvals)

def get_options():
    """
    cosmocalc.py [options] redshift [name_unit [name_unit2 ...]]

    Allowed ``name_unit`` values::

      DA DA_Gyr DA_Mpc DA_cm
      DL DL_Gyr DL_Mpc DL_cm
      DCMR DCMR_Gyr DCMR_Mpc DCMR_cm
      PS_kpc PS_cm
      DTT DTT_Gyr
      VCM VCM_Gpc3
      age age_Gyr
      zage zage_Gyr
      H0 WM WV WK WR z

    If no ``name_unit`` values are supplied then all the above will be printed."""
    from optparse import OptionParser

    parser = OptionParser(get_options.__doc__)
    parser.set_defaults()
    parser.add_option("--H0",
                      default=None,
                      type='float',
                      help="Hubble constant")
    parser.add_option("--WM",
                      default=None,
                      type='float',
                      help="")
    parser.add_option("--WV",
                      default=None,
                      type='float',
                      help="")
    opt, args = parser.parse_args()
    return opt, args, parser

def main():
    opt, args, parser = get_options()

    if len(args) < 1:
        parser.error('Need a redshift')

    kwargs = dict((key, val) for (key, val) in opt.__dict__.items() if val is not None)
    z = float(args[0])
    cc = cosmocalc(z, **kwargs)
    try:
        outlines = []
        for outkey in (args[1:] or _outvals):
            outlines.append(outkey + ' = ' + str(cc[outkey]))
        print('\n'.join(outlines))
    except KeyError:
        parser.error(outkey + ' is not a valid output name_unit')

if __name__ == '__main__':
    main()
