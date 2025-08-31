#ifndef FD_H
#define FD_H

#include <stdio.h>

#include "par.h"

float* ricker(int nt, float dt, float fmax);

void 
fd
( fdFields   *fld, 
  modelPar   *mdl,
  waveletPar *wav,
  geomPar    *geom,
  snapshots  *snap );

#endif // FD_H
