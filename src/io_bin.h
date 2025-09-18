#ifndef BIN_H
#define BIN_H

#include <stdio.h>

void read_receivers(const char *PATH, int f_lines, int f_cols, float *rIdx, float *rIdz);

void read_sources(const char *PATH, int f_lines, int f_cols, float *sIdx, float *sIdz);

void read2D(const char* PATH, void* arr, size_t type, int row, int column);

void write2D(const char* PATH, void* arr, size_t type, int row, int column);

#endif /* BIN_H */
