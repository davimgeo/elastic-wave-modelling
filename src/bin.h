#ifndef BIN_H
#define BIN_H

#include <stdlib.h>

void read_comma_separed_file(const char *path, int nx);

float* read_f32_bin_model(const char *path, int nx, int nz);

void read2D(const char* PATH, void* arr, size_t type, int row, int column);

void write2D(const char* PATH, void* arr, size_t type, int row, int column);

#endif /* BIN_H */
