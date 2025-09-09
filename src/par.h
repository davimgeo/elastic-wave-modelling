#ifndef PAR_H
#define PAR_H

#include <stdio.h>
#include <stdbool.h>

typedef struct 
{
  char *vp_path, *vs_path, *rho_path;

  int   nx, nz, nb;
  float factor;
  float dx, dz;

  float *vp, *vs, *rho;
} modelPar;

typedef struct 
{
  bool snap_bool;
  int  snap_num;
} snapshots;

typedef struct 
{
  float *txx, *tzz, *txz;
  float *vx, *vz;
  
  /* result p wavefield calculated by
   * normalizing txx and tzz */
  float *calc_p;
} fdFields;

typedef struct 
{
  int sIdx, sIdz;
} geomPar;

typedef struct 
{
  int    nt;
  float  dt;
  float  fmax;
  float *wavelet;
} waveletPar;

#endif /* PAR_H */

