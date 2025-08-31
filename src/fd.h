#ifndef FD_H
#define FD_H

#include <stdio.h>

#include "par.h"

float* ricker(int nt, float dt, float fmax);

void 
fd
( fdFields  *fld, 
  float     *vp,
  float     *vs,
  float     *rho,
  int        nx,
  int        nz,
  int        nt,
  float     *wavelet,
  float      dt,
  float      dx,
  float      dz,
  int        sIdx,
  int        sIdz,
  snapshots *snap );

#endif // FD_H
