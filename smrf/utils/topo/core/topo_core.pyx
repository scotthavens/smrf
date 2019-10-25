# cython: embedsignature=True
"""
C implementation of some radiation functions
"""


import cython
import numpy as np
cimport numpy as np
import ctypes


# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()


cdef extern from "topo_core.h":
    void hor1f(int n, double *z, int *h);
    void horval(int n, double *z, double delta, int *h, double *hcos);

@cython.boundscheck(False)
@cython.wraparound(False)
# https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC
def c_hor1d(np.ndarray[double, mode="c", ndim=1] z,
           double spacing,
           np.ndarray[double, mode="c", ndim=1] hcos):
    """
    Call the function hor1f in hor1f.c

    Couldn't figure out how to cast a np.int to a C int so
    I'm just calling it a double

    https://stackoverflow.com/questions/23435756/passing-numpy-integer-array-to-c-code

    Args:
        z: elevation array
        spacing: grid spacing
    
    Returns
        hcos: cosine angle of horizon array changed in place
    """

    cdef int n
    #cdef int m
    n = z.shape[0]
    #m = z.shape[1]

    # convert the z array to C
    cdef np.ndarray[double, mode="c", ndim=1] z_arr
    z_arr = np.ascontiguousarray(z, dtype=np.float64)

    # integer array for horizon index
    cdef np.ndarray[int, ndim=1, mode='c'] h = np.empty((n,), dtype = ctypes.c_int)


    # for i in range(nrows):

    # call the hor1f C function
    hor1f(n, &z_arr[0], &h[0])
    
    # call the horval C function
    horval(n, &z_arr[0], spacing, &h[0], &hcos[0])
