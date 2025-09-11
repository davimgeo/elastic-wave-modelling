#ifndef BIN_H
#define BIN_H

void read_comma_separed_file(const char *path, int nx);

float* read_f32_bin_model(const char *path, int nx, int nz);

void write_f32_bin_model(const char* path, float* model, int nx, int nz);

#endif /* BIN_H */
