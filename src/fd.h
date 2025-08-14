#ifndef FD_H
#define FD_H

#include <stdio.h>

#include "par.h"

float* fd(fdFields* fld, modelPar* model, geomPar* geom, snapshots* snap, waveletPar* wav);

#endif // FD_H