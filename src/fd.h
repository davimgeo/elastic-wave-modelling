#ifndef FD_H
#define FD_H

#include "par.h"

typedef struct 
{
  float *x, *z;
} damping_t;

void fd(const config_t *p, model_t *m, fields_t *fld);

#endif // FD_H
