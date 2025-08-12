#include <stdio.h>
#include <stdlib.h>

#include "read_bin.h"

float* read_f32_bin_model(const char* path, int nx, int nz)
{
    FILE* model_file = fopen(path, "rb"); 
    if (model_file == NULL) {
        printf("Model file could not be opened\n");
        return NULL;
    }

    float *model = (float*)malloc((nx * nz) * sizeof(float));
    for (int i = 0; i < nx; i++) {
        for (int j = 0 ; j < nz; j++) {
            fread(&model[i + j*nz], sizeof(float), 1, model_file); 
        }
    } 

    fclose(model_file);  
    return model;
}

// TODO: correct this function
void write_f32_bin_model(const char* path, float* arr)
{
    FILE* fp = fopen(path, "wb");

    if (fp == NULL) {
        perror("Error oppening file");
        return;
    }

    fwrite(&arr, sizeof(arr), 1, fp);
    fclose(fp);
}
