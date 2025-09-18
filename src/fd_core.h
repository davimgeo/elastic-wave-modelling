#ifndef FD_H
#define FD_H

#include "cfg_par.h"

typedef struct 
{ 
  float *txx;
  float *tzz;
  float *txz;
  float *vx;
  float *vz;
  float *calc_p;
} fields_t;

typedef struct 
{
  float *x, *z;
} damping_t;

void fd(const config_t *cfg, fields_t *fld, model_t *m);

#endif // FD_H
