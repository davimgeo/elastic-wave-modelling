#ifndef READ_BIN_H
#define READ_BIN_H

typedef struct 
{
    int nx;
    int nz;
} m_dim;

float* read_f32_bin_model(const char *path, int nx, int nz);

#endif /* INCLUDE READ_BIN_H */
