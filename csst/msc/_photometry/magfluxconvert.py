# author zouhu
import numpy as np


def fluxerr2magerr(flux, fluxerr, asinh=False, filter='u', zp=22.5):
    """
    convert flux and flux error to mag and mag error (in pogson or asinh form)

    Parameters:
    flux, fluxerr in nanamaggie

    return mag and magerr
    """
    flux = np.array(flux)
    fluxerr = np.array(fluxerr)
    # f0=1.0e9
    f0 = 10 ** (zp / 2.5)
    nn = flux.size
    mag = np.array(np.ones_like(flux) * 99.0)
    magerr = np.array(np.ones_like(flux) * 99.0)
    if not asinh:
        mask = flux > 0
        if mask.any():
            mag[mask] = -2.5 * np.log10(flux[mask] / f0)
            magerr[mask] = 2.5 / np.log(10.0) * fluxerr[mask] / flux[mask]
    else:
        bs = {'u': 1.4e-10, 'g': 0.9e-10, 'r': 1.2e-10, 'i': 1.8e-10, 'z': 7.4e-10}
        b = bs[filter]
        mag = -(2.5 / np.log(10.0)) * (np.arcsinh((flux / f0) / (2.0 * b)) + np.log(b))
        magerr = 2.5 / np.log(10.0) * (fluxerr / f0) / (2.0 * b) / np.sqrt(1.0 + ((flux / f0) / (2.0 * b)) ** 2)
    return mag, magerr


def magerr2fluxerr(mag, magerr, asinh=False, filter='u', zp=22.5):
    """
    convert mag and mag error to flux and flux error (in pogson or asinh form)

    Parameters:
        mag,magerr ndarray

    return
        flux, fluxerr in nanamaggie
    """
    # f0=1.0e9
    f0 = 10 ** (zp / 2.5)
    if not asinh:
        flux = 10.0 ** (mag / (-2.5)) * f0
        fluxerr = flux * magerr * np.log(10.0) / 2.5
    else:
        bs = {'u': 1.4e-10, 'g': 0.9e-10, 'r': 1.2e-10, 'i': 1.8e-10, 'z': 7.4e-10}
        b = bs[filter]
        flux = np.sinh(mag / (-2.5 / np.log(10.0)) - np.log(b)) * 2.0 * b * f0
        fluxerr = magerr * np.log(10.0) / 2.5 * (2.0 * b) * np.sqrt(1.0 + ((flux / f0) / (2.0 * b)) ** 2) * f0
    return flux, fluxerr


def asinhpogson(mag1, magerr1, asinh2pogson=False, filter='u', zp=22.5):
    """
    convert magnitude form between asinh and pogson
    """
    flux, fluxerr = magerr2fluxerr(mag1, magerr1, asinh=asinh2pogson, filter=filter, zp=zp)
    mag2, magerr2 = fluxerr2magerr(flux, fluxerr, asinh=(not asinh2pogson), filter=filter, zp=zp)
    return mag2, magerr2
