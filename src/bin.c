#include <stdio.h>
#include <stdlib.h>

#include "bin.h"

float* read_f32_bin_model(const char* path, int nx, int nz)
{
    FILE* model_file = fopen(path, "rb"); 
    if (model_file == NULL) {
        printf("Model file could not be opened\n");
        return NULL;
    }

    size_t n = nx * nz;
    float *model = (float*)calloc(n, sizeof(float));
    for (int i = 0; i < nx; i++) {
        for (int j = 0 ; j < nz; j++) {
            fread(&model[i + j*nz], sizeof(float), 1, model_file); 
        }
    } 

    fclose(model_file);  
    return model;
}

void write_f32_bin_model(const char* path, float* model, int nx, int nz)
{
    FILE* model_file = fopen(path, "wb");
    if (model_file == NULL) {
        printf("Model file could not be created\n");
        return;
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < nz; j++) {
            fwrite(&model[i + j*nz], sizeof(float), 1, model_file);
        }
    }

    fclose(model_file);
}
