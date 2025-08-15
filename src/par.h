#ifndef PAR_H
#define PAR_H

#include <stdio.h>
#include <stdbool.h>

typedef struct {
  const char *rho_path;
  const char *vs_path;
  const char *vp_path;
  int         nx;
  int         nz;
  float       dx;
  float       dz;

  float *     rho;
  float *     vs;
  float *     vp;
} modelPar;

typedef struct {
  bool snap_bool;
  int  snap_num;
} snapshots;

typedef struct {
  float *txx;
  float *tzz;
  float *txz;
  float *vx;
  float *vz;
} fdFields;

typedef struct {
  int sIdx;
  int sIdz;
} geomPar;

typedef struct {
  int    nt;
  float  dt;
  float  fmax;
  float *wavelet;
} waveletPar;

#endif /* PAR_H */

