import numpy as np
from scipy import ndimage

from smrf.utils.topo import hor1d

import matplotlib.pyplot as plt


def viewf(dem, spacing, nangles=16, slope=None):
    """
    Calculate the sky view factor of a dem

    Args:
        spacing:
        dem:
        nangles:

    """

    if dem.ndim != 2:
        raise ValueError('viewf input of dem is not a 2D array')

    if nangles != 16 and nangles != 32:
        raise ValueError('viewf number of angles can be 16 or 32')

    hcos = {}

    # East
    hcos['e'] = hor1d.hor2d_c(dem, spacing)

    # West
    hcos['w'] = hor1d.hor2d_c(np.ascontiguousarray(dem[::-1]), spacing)

    # SSW
    t = skew(dem, -22.5).transpose()
    hcos['ssw'] = skew(hor1d.hor2d_c(
        np.ascontiguousarray(t), spacing).transpose(), -22.5, fwd=False)

    # NNE
    hcos['nne'] = skew(hor1d.hor2d_c(
        np.ascontiguousarray(t[::-1]), spacing).transpose(), -22.5, fwd=False)

    # SW
    t = skew(dem, -45).transpose()
    hcos['sw'] = skew(hor1d.hor2d_c(np.ascontiguousarray(t),
                                    spacing).transpose(), -45, fwd=False)

    # NE
    hcos['ne'] = skew(hor1d.hor2d_c(np.ascontiguousarray(
        t[::-1]), spacing).transpose(), -45, fwd=False)

    # SSE
    t = skew(dem, 22.5).transpose()
    hcos['sse'] = skew(hor1d.hor2d_c(np.ascontiguousarray(t),
                                     spacing).transpose(), 22.5, fwd=False)

    # NNW
    hcos['nnw'] = skew(hor1d.hor2d_c(np.ascontiguousarray(
        t[::-1]), spacing).transpose(), 22.5, fwd=False)

    # SE
    t = skew(dem, 45).transpose()
    hcos['se'] = skew(hor1d.hor2d_c(np.ascontiguousarray(t),
                                    spacing).transpose(), 45, fwd=False)

    # NW
    hcos['nw'] = skew(hor1d.hor2d_c(np.ascontiguousarray(
        t[::-1]), spacing).transpose(), 45, fwd=False)

    # S
    demt = dem.transpose()
    hcos['s'] = hor1d.hor2d_c(np.ascontiguousarray(demt),
                              spacing).transpose()

    # N
    hcos['n'] = hor1d.hor2d_c(np.ascontiguousarray(demt[::1]),
                              spacing).transpose()

    # ENE
    t = skew(demt, -22.5).transpose()
    hcos['ene'] = skew(hor1d.hor2d_c(np.ascontiguousarray(t),
                                     spacing).transpose(), -22.5, fwd=False).transpose()

    # WSW
    hcos['wsw'] = skew(hor1d.hor2d_c(np.ascontiguousarray(
        demt[::-1]), spacing).transpose(), -22.5, fwd=False).transpose()

    # ESE
    t = skew(demt, 22.5).transpose()
    hcos['ese'] = skew(hor1d.hor2d_c(np.ascontiguousarray(t),
                                     spacing).transpose(), 22.5, fwd=False).transpose()

    # WNW
    hcos['wnw'] = skew(hor1d.hor2d_c(np.ascontiguousarray(
        demt[::-1]), spacing).transpose(), 22.5, fwd=False).transpose()

    # calcualte the gradient
    if slope is None:
        pass

    svf = viewcalc(slope, hcos)

    return svf


def viewcalc(gradient, hcos):
    """
    """
    pass


def skew(arr, angle, fwd=True):
    """
    Skew the origin of successive lines by a specified angle
    A skew with angle of 30 degrees causes the following transformation:

        +-----------+       +---------------+
        |           |       |000/          /|
        |   input   |       |00/  output  /0|
        |   image   |       |0/   image  /00|
        |           |       |/          /000|
        +-----------+       +---------------+

    Calling skew with fwd=False will return the output image
    back to the input image.

    Skew angle must be between -45 and 45 degrees

    Args:
        arr: array to skew
        angle: angle between -45 and 45 to skew by
        fwd: add skew to image if True, unskew image if False

    Returns:
        skewed array

    """

    if angle == 0:
        return arr

    if angle > 45 or angle < -45:
        raise ValueError('skew angle must be between -45 and 45 degrees')

    nlines, nsamps = arr.shape

    if angle >= 0.0:
        negflag = False
    else:
        negflag = True
        angle = -angle

    slope = np.tan(angle * np.pi / 180.0)
    max_skew = int((nlines - 1) * slope + 0.5)

    o_nsamps = nsamps
    if fwd:
        o_nsamps += max_skew
    else:
        o_nsamps -= max_skew

    b = np.zeros((nlines, o_nsamps))
    for line in range(nlines):
        o = line if negflag else nlines - line - 1
        offset = int(o * slope + 0.5)

        if fwd:
            b[line, offset:offset+nsamps] = arr[line, :]
        else:
            b[line, :] = arr[line, offset:offset+o_nsamps]

    return b
