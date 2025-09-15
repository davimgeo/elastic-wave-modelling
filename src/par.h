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

  /* snapshot parameters */
  bool snap_bool;
  int  snap_num;
  int  snap_ratio;

  /* geometry parameters */
  int sIdx, sIdz;

  /* wavelet parameters */
  int    nt;
  float  dt;
  float  fmax;
  float *wavelet;

} config_t;

typedef struct
{
  float *vp, *vs, *rho;
} model_t;


typedef struct 
{ 
  float *txx, *tzz, *txz;
  float *vx, *vz;
  float *calc_p;
} fields_t;

#endif /* PAR_H */

