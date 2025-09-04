#ifndef FD_H
#define FD_H

#include "par.h"

float* ricker(int nt, float dt, float fmax);

void set_boundary(fdFields *fld, modelPar *mld);

void 
fd
( fdFields   *fld, 
  modelPar   *mdl,
  waveletPar *wav,
  geomPar    *geom,
  snapshots  *snap );

#endif // FD_H
