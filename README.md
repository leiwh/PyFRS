## PyFS
This python code is dedicated to compute the multiband synchrotron emissions from the Forward-Reverse Shocks of a GRB or TDE jet.


## Features
fs.py computes synchrotron emission from the forward shock of a jet/outflow. It includes:

 - Synchrotron emission from relativistic jet or un-relativistic outflow;
 - Jet/outflow dynamics cover Coasting, deceleration and Newtonian phases;
 - Arbitray medium density profile;
 - Arbitray viewing angles;
 - Self-absorbtion, SSC and EIC;
 - Jet break effect;
 - Energy injection;
 - Density jump.

rs.py computes synchrotron emission from the reverse shock of a jet/outflow.

Future:

 - Neutrino emission

## Using
In your python code, import the library with `import fs`.  

The main function of interest is`fs.FS_flux(t, nu, **Z)`.  See `LightCurve.ipynb` for a simple example.

1. t: a 1-D array of observer-frame times, in seconds (s)
2. nu: a 1-D array, the same shape as t, of observer-frame frequencies, in Hertz (Hz)

3. Fnu is an array, same size as t and nu, containing the observed flux in mJy at each time and frequency.

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



# References
Lei, W.-H., Yuan, Q., Zhang, B., & Wang, D. 2016, ApJ, 816, 20