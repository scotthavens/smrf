import numpy as np
from math import hypot

from smrf.utils.topo.core import topo_core


def hor2d_c(z, spacing, fwd=True):
    """
    Calculate values of cosines of angles to horizons in 2 dimension, 
    measured from zenith, from elevation difference and distance.  Let
    G be the horizon angle from horizontal and note that:

        sin G = z / sqrt( z^2 + dis^2);

    This result is the same as cos H, where H measured from zenith.

    Args:
        z: elevation array
        spacing: spacing of array

    Returns:
        hcos: cosines of angles to horizon
    """

    if z.ndim != 2:
        raise ValueError('hor1d input of z is not a 2D array')

    if z.dtype != np.double:
        raise ValueError('hor1d input of z must be a double')

    spacing = np.double(spacing)

    if not fwd:
        z = np.ascontiguousarray(z[::-1])
    else:
        z = np.ascontiguousarray(z)

    h = np.zeros_like(z)
    topo_core.c_hor2d(z, spacing, h)

    if not fwd:
        h = h[::-1]

    return h


def hor1d_c(z, spacing):
    """
    Calculate values of cosines of angles to horizons in 1 dimension, 
    measured from zenith, from elevation difference and distance.  Let
    G be the horizon angle from horizontal and note that:

        sin G = z / sqrt( z^2 + dis^2);

    This result is the same as cos H, where H measured from zenith.

    Args:
        z: elevation array
        spacing: spacing of array

    Returns:
        hcos: cosines of angles to horizon
    """

    if z.ndim != 1:
        raise ValueError('hor1d input of z is not a 1D array')

    h = np.zeros_like(z)
    topo_core.c_hor1d(z, spacing, h)

    return h


def hor1d(z, spacing):
    """
    Calculate values of cosines of angles to horizons, measured
    from zenith, from elevation difference and distance.  Let
    G be the horizon angle from horizontal and note that:

        sin G = z / sqrt( z^2 + dis^2);

    This result is the same as cos H, where H measured from zenith.

    Args:
        z: elevation array
        spacing: spacing of array

    Returns:
        hcos: cosines of angles to horizon
    """

    if z.ndim != 2:
        raise ValueError('hor1d input of z is not a 2D array')

    hcos = np.zeros_like(z)
    nrow = z.shape[0]

    for i in range(nrow):
        h = hor1f(z[i, :])
        hcos[i, :] = horval(z[i, :], spacing, h)

    return hcos


def horval(z, delta, h):
    """
    Calculate values of cosines of angles to horizons, measured
    from zenith, from elevation difference and distance.  Let
    G be the horizon angle from horizontal and note that:

        sin G = z / sqrt( z^2 + dis^2);

    This result is the same as cos H, where H measured from zenith.

    Args:
        z: elevation vector
        delta: spacing
        h: horizon function output

    Returns:
        hcos: cosines of angles to horizon
    """

    hcos = np.zeros_like(z)
    for i in range(len(z)):

        # grid points to horizon
        j = h[i]
        d = j - i

        # point is its own horizon
        if (d == 0):
            hcos[i] = 0

        # else need to calculate sine */
        else:
            if (d < 0):
                d = -d

            diff = z[j] - z[i]
            hcos[i] = diff / hypot(diff, d * delta)

    return hcos


def hor1f(z):
    """
    Calculate the horizon pixel for all z
    This mimics the algorthim from Dozier 1981 and the
    hor1f.c from IPW

    Works backwards from the end but looks forwards for
    the horizon

    xrange stops one index before [stop]

    Args:
        z - elevations for the points

    Returns:
        h - index to the horizon point

    20150601 Scott Havens
    20191025 Scott Havens
    """

    N = len(z)  # number of points to look at
    # x = np.array(x)
    z = np.array(z)

    # preallocate the h array to zeros, what ealloc is doing
    h = np.zeros(N, dtype=int)

    # the end point is it's own horizon
    h[N-1] = N-1

    # loop runs from next-to-end backwards to the beginning
    # range end is -1 to get the 0 index
    for i in range(N-2, -1, -1):

        zi = z[i]

        # Start with next-to-adjacent point in either forward or backward
        # direction, depending on which way loop is running. Note that we
        # don't consider the adjacent point; this seems to help reduce noise.
        k = i + 2

        if k >= N:
            k -= 1

        # loop until horizon is found
        # xrange will set the maximum number of iterations that can be
        # performed based on the length of the vector

        j = k
        k = h[j]
        sij = _slope(i, zi, j, z[j])
        sihj = _slope(i, zi, k, z[k])

        # if slope(i,j) >= slope(i,h[j]), horizon has been found; otherwise
        # set j to k (=h[j]) and loop again
        # or if we are at the end of the section
        while sij < sihj:

            j = k
            k = h[j]

            sij = _slope(i, zi, j, z[j])
            sihj = _slope(i, zi, k, z[k])

        # if slope(i,j) > slope(j,h[j]), j is i's horizon; else if slope(i,j)
        # is zero, i is its own horizon; otherwise slope(i,j) = slope(i,h[j])
        # so h[j] is i's horizon
        if sij > sihj:
            h[i] = j
        elif sij == 0:
            h[i] = i
        else:
            h[i] = k

    return h


def hor1f_vec(z):
    """
    Calculate the horizon pixel for all z
    This trys to improve the method by vectorizing it
    with Numpy as for loops are slow

    Works backwards from the end but looks forwards for
    the horizon

    Args:
        z - elevations for the points

    Returns:
        h - index to the horizon point

    20150601 Scott Havens
    20191025 Scott Havens
    """

    N = len(z)  # number of points to look at
    # x = np.array(x)
    z = np.array(z)

    # preallocate the h array to zeros, what ealloc is doing
    h = np.zeros(N, dtype=int)

    # the end point is it's own horizon
    h[N-1] = N-1

    # loop runs from next-to-end backwards to the beginning
    # range end is -1 to get the 0 index
    for i in range(N-2, -1, -1):

        zi = z[i]

        # Start with next-to-adjacent point in either forward or backward
        # direction, depending on which way loop is running. Note that we
        # don't consider the adjacent point; this seems to help reduce noise.
        k = i + 2

        if k >= N:
            k -= 1

        # loop until horizon is found
        # xrange will set the maximum number of iterations that can be
        # performed based on the length of the vector

        j = k
        k = h[j]
        sij = _slope(i, zi, j, z[j])
        sihj = _slope(i, zi, k, z[k])

        # if slope(i,j) >= slope(i,h[j]), horizon has been found; otherwise
        # set j to k (=h[j]) and loop again
        # or if we are at the end of the section
        while sij < sihj:

            j = k
            k = h[j]

            sij = _slope(i, zi, j, z[j])
            sihj = _slope(i, zi, k, z[k])

        # if slope(i,j) > slope(j,h[j]), j is i's horizon; else if slope(i,j)
        # is zero, i is its own horizon; otherwise slope(i,j) = slope(i,h[j])
        # so h[j] is i's horizon
        if sij > sihj:
            h[i] = j
        elif sij == 0:
            h[i] = i
        else:
            h[i] = k

    return h


def _slope(xi, zi, xj, zj):
    """
    Slope between the two points only if the pixel is higher
    than the other
    20150603 Scott Havens
    """

    if zj <= zi:
        return 0
    else:
        return (zj - zi) / (xj - xi)
