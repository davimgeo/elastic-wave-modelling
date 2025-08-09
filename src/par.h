#ifndef PAR_H 
#define PAR_H

#include <stdio.h>

typedef struct _modelPar {
    const char* path;
    int nx;
    int nz;
} modelPar;

typedef struct _waveletPar {
    int nt;
    float dt;
    float fmax;
} waveletPar;

#endif // PAR_H