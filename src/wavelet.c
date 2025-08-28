#include <math.h>
#include <stdlib.h>

#define PI 3.1415926f

float* ricker(int nt, float dt, float fmax) 
{
  float* ricker = (float *)calloc((size_t)nt, sizeof(float));
  if (ricker == NULL) {
    return NULL;
  }
  
  float t0 = 2.0f * PI / fmax;
  float fc = fmax / (3.0f*sqrtf(PI));

  for (int i = 0; i < nt; i++) {
    float t = (i * dt) - t0;
    float arg = PI*(PI*PI * fc*fc * t*t);
    ricker[i] = (1.0f - 2.0f*arg) * expf(-arg);
  }
  return ricker;
}

