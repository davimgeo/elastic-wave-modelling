#ifndef BIN_H
#define BIN_H

#include <stdlib.h>

void read2D(const char* PATH, void* arr, size_t type, int row, int column);

void write2D(const char* PATH, void* arr, size_t type, int row, int column);

#endif /* BIN_H */
