from numpy cimport ndarray
from numpy import empty

cimport cdirac

def rdirac(double sol, int n, int l, int k, int np,
        ndarray[double, mode="c"] r not None,
        ndarray[double, mode="c"] vr not None,
        double E):
    cdef int nr = len(r)
    assert len(vr) == nr
    cdef ndarray[double, mode="c"] R = empty(nr, dtype="double")
    cdirac.c_rdirac(&sol, &n, &l, &k, &np, &nr, &r[0], &vr[0], &E,
            &R[0])
    return E, R
