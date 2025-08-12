#ifndef PAR_H 
#define PAR_H

#include <stdio.h>

typedef struct _modelPar {
    const char* path;
    int   nx;
    int   nz;
    float dx;
    float dz;
} modelPar;

typedef struct _geomPar {
    int srcX;
} geomPar;

typedef struct _waveletPar {
    int nt;
    float dt;
    float fmax;
} waveletPar;

#endif // PAR_H
