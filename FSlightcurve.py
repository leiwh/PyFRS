import numpy as np
import matplotlib.pyplot as plt
try:
    from PyFRS import cgs,fs,rs
except ImportError:
    import cgs,fs,rs



def FS_lc(v,**Z):

    t=np.geomspace(1.e1,1.e8,500)
    nu=np.empty(t.shape)
    nu[:]=v

    Fnu=fs.FS_flux(t,nu,**Z)

    fig,ax=plt.subplots(1,1)
    ax.plot(t,Fnu)
    ax.set(xscale='log', xlabel=r'$t$ (s)',
           yscale='log', ylabel=r'$F_\nu$[$10^{9}$ Hz] (mJy)')

    fig.tight_layout()
    print("Saving figure lc.png")
    fig.savefig("FS_lc.png")
    plt.close(fig)
