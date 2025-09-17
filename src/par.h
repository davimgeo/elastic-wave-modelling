#ifndef PAR_H
#define PAR_H

typedef struct
{
  char  *vp_path, *vs_path, *rho_path;

  int    nx, nz, nb;
  float  dx, dz;

  float  factor;
  int    nxx, nzz;
} model_par_t;

typedef struct 
{
  float *vp;
  float *vs;
  float *rho;
} model_t;

typedef struct
{
  int snap_bool;
  int snap_num;
  int snap_ratio;
} snap_t;

typedef struct
{
  float *sIdx;
  float *sIdz;
  float *rIdx; 
  float *rIdz;
} geom_t;

typedef struct
{
  int    nt;
  float  dt;
  float  fmax;
  float *wavelet;
} wavelet_t;

#endif // PAR_H
