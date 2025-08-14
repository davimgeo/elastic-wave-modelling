#ifndef BIN_H
#define BIN_H

float* read_f32_bin_model(const char *path, int nx, int nz);

void write_f32_bin_model(const char* path, float* model, int nx, int nz);

#endif /* BIN_H */
