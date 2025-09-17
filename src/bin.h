#ifndef BIN_H
#define BIN_H

#include <stdlib.h>

#include "par.h"

void read_geometry(const char *path, int f_lines, int f_cols, geom_t *geom);

void read2D(const char* PATH, void* arr, size_t type, int row, int column);

void write2D(const char* PATH, void* arr, size_t type, int row, int column);

#endif /* BIN_H */
