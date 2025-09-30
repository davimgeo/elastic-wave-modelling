#ifndef FD_H
#define FD_H

#include "par.h"

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

void fd_normal(const config_t *p, model_t *m, fields_t *fld);

void fd_simd(const config_t *p, model_t *m, fields_t *fld);

#endif // FD_H
