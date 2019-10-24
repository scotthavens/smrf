import numpy as np


def ihorizon(x, y, Z, azm, mu=0, offset=2, ncores=0):
    """
    Calculate the horizon values for an entire DEM image
    for the desired azimuth

    Assumes that the step size is constant

    Args:
        X - vector of x-coordinates
        Y - vector of y-coordinates
        Z - matrix of elevation data
        azm - azimuth to calculate the horizon at
        mu - 0 -> calculate cos(z)
             - >0 -> calculate a mask whether or not the point can see the sun

    Returns:
        H   - if mask=0 cosine of the local horizonal angles
            - if mask=1 index along line to the point

    20150602 Scott Havens
    """

    # check inputs
    azm = azm*np.pi/180  # degress to radians
    m, n = Z.shape

    # transform the x,y into the azm direction xr,yr
    xi, yi = np.arange(-n/2, n/2), np.arange(-m/2, m/2)
    X, Y = np.meshgrid(xi, yi)
    xr = X*np.cos(azm) - Y*np.sin(azm)
    yr = X*np.sin(azm) + Y*np.cos(azm)

    # xr is the "new" column index for the profiles
    # yr is the distance along the profile
    xr = xr.round().astype(int)
    yr = (x[2] - x[1]) * yr

    H = np.zeros(Z.shape)
#     pbar = progressbar.ProgressBar(n).start()
#     j = 0

    # loop through the columns
#     if ncores == 0:
    for i in np.arange(-n/2, n/2):
        find_horizon(i, H, xr, yr, Z, mu)

#     else:
#     shared_array_base = mp.Array(ctypes.c_double, m*n)
#     sH = np.ctypeslib.as_array(shared_array_base.get_obj())
#     sH = sH.reshape(m, n)
#     sxr = np.ctypeslib.as_array(shared_array_base.get_obj())
#     sxr = sxr.reshape(m, n)
#     syr = np.ctypeslib.as_array(shared_array_base.get_obj())
#     syr = syr.reshape(m, n)
#     sZ = np.ctypeslib.as_array(shared_array_base.get_obj())
#     sZ = sZ.reshape(m, n)
#     def wrap_horizon(i, def_param1=sH, def_parm2=sxr, def_param3=syr, def_param4=sZ):
#         find_horizon(i, sH, sxr, syr, sZ, 0.67)

#     pool = mp.Pool(processes=4)
#     [pool.apply(find_horizon, args=(i, shared_array, xr, yr, Z, mu, offset)) for i in range(-n/2,n/2)]
#     pool.map(wrap_horizon, range(-n/2,n/2))

#     print(shared_array)

    return H


def find_horizon(i, H, xr, yr, Z, mu):
    # index to profile and get the elevations
    ind = xr == i
    zi = Z[ind]

    # distance along the profile
    di = yr[ind]

    # sort the y values and get the cooresponding elevation
    idx = np.argsort(di)
    di = di[idx]
    zi = zi[idx]

    # if there are some values in the vector
    # calculate the horizons
    if len(zi) > 0:
        # h2 = hor1f_simple(di, zi)
        h = hor1f(di, zi)

        cz = _cosz(di, zi, di[h], zi[h])

        # if we are making a mask
        if mu > 0:
            # iz = cz == 0    # points that are their own horizon
            idx = cz > mu   # points sheltered from the sun
            cz[idx] = 0
            cz[~idx] = 1
#                 cz[iz] = 1

    H[ind] = cz

#         j += 1
#         pbar.update(j)
#     pbar.finish()


def hord(z):
    """
    Calculate the horizon pixel for all x,z
    This mimics the simple algorthim from Dozier 1981
    to help understand how it's working

    Works backwards from the end but looks forwards for
    the horizon 90% faster than rad.horizon

    Args::
        x - horizontal distances for points
        z - elevations for the points

    Returns:
        h - index to the horizon point

    20150601 Scott Havens
    """

    N = len(z)  # number of points to look at
#     offset = 1      # offset from current point to start looking

    # preallocate the h array
    h = np.zeros(N, dtype=int)
    h[N-1] = N-1
    i = N - 2

    # work backwarks from the end for the pixels
    while i >= 0:
        h[i] = i
        j = i + 1   # looking forward
        max_tan = -9999

        while j < N:
            sij = _slope_all(i, z[i], j, z[j])

            if sij > max_tan:
                h[i] = j
                max_tan = sij

            j = j + 1
        i = i - 1

    return h


def hor1f_simple(z):
    """
    Calculate the horizon pixel for all x,z
    This mimics the simple algorthim from Dozier 1981
    to help understand how it's working

    Works backwards from the end but looks forwards for
    the horizon 90% faster than rad.horizon

    Args:
        x - horizontal distances for points
        z - elevations for the points

    Returns:
        h - index to the horizon point

    20150601 Scott Havens
    """

    N = len(z)  # number of points to look at
#     offset = 1      # offset from current point to start looking

    # preallocate the h array
    h = np.zeros(N, dtype=int)
    h[N-1] = N-1
    i = N - 2

    # work backwarks from the end for the pixels
    while i >= 0:
        h[i] = i
        j = i + 1   # looking forward
        max_tan = 0

        while j < N:
            sij = _slope(i, z[i], j, z[j])

            if sij > max_tan:
                h[i] = j
                max_tan = sij

            j = j + 1
        i = i - 1

    return h


def hor1f(x, z, offset=1):
    """
    BROKEN: Haven't quite figured this one out

    Calculate the horizon pixel for all x,z
    This mimics the algorthim from Dozier 1981 and the
    hor1f.c from IPW

    Works backwards from the end but looks forwards for
    the horizon

    xrange stops one index before [stop]

    Args:
        x - horizontal distances for points
        z - elevations for the points

    Returns:
        h - index to the horizon point

    20150601 Scott Havens
    """

    N = len(x)  # number of points to look at
    x = np.array(x)
    z = np.array(z)

    # preallocate the h array
    h = np.zeros(N, dtype=int)
    h[N-1] = N-1    # the end point is it's own horizon

    # work backwards from the end for the pixels
    for i in range(N-2, -1, -1):

        zi = z[i]

        # Start with next-to-adjacent point in either forward or backward
        # direction, depending on which way loop is running. Note that we
        # don't consider the adjacent point; this seems to help reduce noise.
        k = i + offset

        if k >= N:
            k -= 1

        # loop until horizon is found
        # xrange will set the maximum number of iterations that can be
        # performed based on the length of the vector
        for t in range(k, N):
            j = k
            k = h[j]

            sij = _slope(x[i], zi, x[j], z[j])
            sihj = _slope(x[i], zi, x[k], z[k])

            # if slope(i,j) >= slope(i,h[j]), horizon has been found; otherwise
            # set j to k (=h[j]) and loop again
            # or if we are at the end of the section
            if sij > sihj:  # or k == N-1:
                break

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

    return 0 if zj <= zi else (zj - zi) / (xj - float(xi))


def _slope_all(xi, zi, xj, zj):
    """
    Slope between the two points only if the pixel is higher
    than the other
    20150603 Scott Havens
    """

    return (zj - zi) / (xj - float(xi))


def _cosz(x1, z1, x2, z2):
    """
    Cosize of the zenith between two points

    20150601 Scott Havens
    """
    d = np.sqrt((x2 - x1)**2 + (z2 - z1)**2)
    diff = z2 - z1

#     v = np.where(diff != 0., d/diff, 100)

    i = d == 0
    d[i] = 1
    v = diff/d
    v[i] = 0

    return v
