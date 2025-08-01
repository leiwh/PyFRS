## PyFRS
This python code is dedicated to compute the multi-band synchrotron emissions from the Forward-Reverse Shocks of a GRB or TDE jet. It is designed based on the book `The Physics of Gamma-Ray Bursts` by Bing Zhang (2018) and the review paper `A complete reference of the analytical synchrotron external shock models of gamma-ray bursts` by Gao et al. (2013, NewAR, 57, 141). The jet dynamics is calculated numerically using the equations given by Huang, Gou, Dai & Lu (2000, ApJ, 543, 90). 

You can use the code in a similar way to afterglowpy (https://github.com/geoffryan/afterglowpy). See Lightcurve.ipynb for an example. The code is lightweight. One can try jetsimpy (https://github.com/haowang-astro/jetsimpy) by Hao Wang for a version based on hydrodynamic simulations, and VegasAfterglow (https://github.com/YihanWangAstro/VegasAfterglow) by Yihan Wang for a more comprehensive one. 


## Features
`fs.py` computes synchrotron emission from the forward shock of a jet/outflow. It includes:

 - Synchrotron emission from relativistic jet or non-relativistic outflow;
 - Jet/outflow dynamics cover Coasting, Deceleration and Newtonian phases;
 - Arbitrary medium density profile;
 - Arbitrary viewing angles (on and off-beam);
 - Self-absorbtion, SSC and EIC;
 - Jet break effect;
 - Energy injection;
 - Density jump.

`rs.py` computes synchrotron emission from the reverse shock of a jet/outflow.
- rs.py is updated version, rs_scaling.py is old version based on simple scalings;
- rs_scaling.py includes only thin shell case.

Future:

 - Neutrino emission


## Installing

Clone the repository:
```
git clone https://github.com/leiwh/PyFRS
```

Or install with pip:
```
pip install PyFRS
```

It runs properly with anaconda/edm and VS Code.

## Usage
The forward and reverse shock modules can be used separately.

### Forward Shock (FS)
For forward shock (FS), import the library with `import fs` in your python code.  

The main function of interest is`fs.FS_flux(t, nu, **Z)`.  See `LightCurve.ipynb` for a simple example.

1. t: a 1-D array of observer-frame times, in seconds (s)
2. nu: a 1-D array, the same shape as t, of observer-frame frequencies, in Hertz (Hz)

3. Fnu=fs.FS_flux(t, nu, **Z) is an array, same size as t and nu, the model flux in mJy at each time and frequency.

4. For forward shock afterglows `Z` has 12 required keyword arguments:

- `z`            Redshift
- `k`            Medium density profile index
- `SSC`          SSC
- `XIC`          EIC
- `E0`           Isotropic-equivalent energy in erg
- `Gamma0`       Initial bulk Lorentz factor
- `theta_j`      Jet opening angle
- `theta_obs`    Viewing angle in deg
- `n18`          CNM density at 10^18cm
- `p`            Electron energy distribution index
- `epsilon_e`    Fraction of the shock energy into electrons
- `epsilon_B`    Fraction of the shock energy in magnetic field


### Reverse Shock (RS)
For reverse shock (RS), import the library with `import rs` in your python code.  

The main function of interest is`rs.RS_flux(t, nu, **Z)`.  See `LightCurve.ipynb` for a simple example.

1. t: a 1-D array of observer-frame times, in seconds (s)
2. nu: a 1-D array, the same shape as t, of observer-frame frequencies, in Hertz (Hz)

3. Fnu=rs.RS_flux(t, nu, **Z) is an array, same size as t and nu, the model flux in mJy at each time and frequency.

4. For reverse shock afterglows `Z` has 12 required keyword arguments:

- `z`            Redshift
- `k`            Medium density profile index
- `SSC`          SSC
- `XIC`          EIC
- `E0`           Isotropic-equivalent energy in erg
- `Gamma0`       Initial bulk Lorentz factor
- `theta_j`      Jet opening angle
- `theta_obs`    Viewing angle in deg
- `n18`          CNM density at 10^18cm
- `p`            Electron energy distribution index
- `epsilon_e`    Fraction of the shock energy into electrons
- `epsilon_B`    Fraction of the shock energy in magnetic field


# References

Zhang, B. 2018, The Physics of Gamma-Ray Bursts

Gao, H., Lei, W.-H., Zou, Y.-C., Wu, X.-F., & Zhang, B. 2013, NewAR, 57, 141

Lei, W.-H., Yuan, Q., Zhang, B., & Wang, D. 2016, ApJ, 816, 20

Zhu, Z.-P., Xu, D., Fynbo, J. P. U., et al. 2023, ApJ, 948, 30