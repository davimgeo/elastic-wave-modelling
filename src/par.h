#ifndef PAR_H
#define PAR_H

#include <stdio.h>
#include <stdbool.h>

typedef struct {
    /* model parameters */
    char  *vp_path, *vs_path, *rho_path;

    int    nx, nz, nb;
    float  dx, dz;

    float  factor;
    int    nxx, nzz;

    float *vp, *vs, *rho;

    /* snapshot parameters */
    bool snap_bool;
    int  snap_num;
    int  snap_ratio;

    /* fields */
    float *txx, *tzz, *txz;
    float *vx, *vz;
    // resulted p measured by 0.5*(txx+tzz)
    float *calc_p;

    /* geometry parameters */
    int sIdx, sIdz;

    /* wavelet parameters */
    int    nt;
    float  dt;
    float  fmax;
    float *wavelet;

} config_t;

#endif /* PAR_H */

