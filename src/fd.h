#ifndef FD_H
#define FD_H

#include "par.h"

typedef struct 
{
  float *x, *z;
} damping_t;

float* ricker(int nt, float dt, float fmax);

void fd(config_t *config);

#endif // FD_H
