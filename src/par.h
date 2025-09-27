#ifndef PAR_H
#define PAR_H

#include <stdio.h>

typedef struct
{
  char  *vp_path, *vs_path, *rho_path;
  int    nx, nz, nb;
  float  dx, dz;
  float  factor;
  int    nxx, nzz;

  int snap_bool;
  int snap_num;
  int snap_ratio;

  int src_f_lines;
  int src_f_cols;

  float *src_x;
  float *src_z;
  int r_f_lines;
  int r_f_cols;

  float *rcv_x; 
  float *rcv_z;

  int    nt;
  float  dt;
  float  fmax;
  float *wavelet;

} config_t;

typedef struct 
{
  float *vp;
  float *vs;
  float *rho;
} model_t;

#endif /* PAR_H */

