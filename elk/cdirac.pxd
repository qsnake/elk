cdef extern:
    void c_rdirac(double *sol, int *n, int *l, int *k, int *np, int *nr,
            double *r, double *vr, double *_eval, double *y)
