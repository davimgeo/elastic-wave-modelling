#ifndef FD_H
#define FD_H

#include "par.h"

typedef struct 
{
  float *x; 
  float *z;
} damping_t;

float* ricker(int nt, float dt, float fmax);

void fd(config_t *p);

#endif // FD_H
