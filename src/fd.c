#include <stdlib.h>
#include <math.h>

#include "fd.h"

#define PI 3.14159265358979f

float* ricker(int nt, float dt, float fmax) 
{
    float* ricker = (float*)malloc(nt * sizeof(float));
    if (ricker == NULL) {
        return NULL;
    }
    
    float t0 = 2.0f * PI / fmax;
    for (int i = 0; i < nt; i++) {
        float t = (i * dt) - t0;
        float arg = (PI*PI * fmax*fmax * t*t);
        ricker[i] = (1.0f - 2.0f*arg) * expf(-arg);
    }
    return ricker;
}

/* implicit declaration */
static void fd2d_2E2T(int nx, int nz, float dx, float dz, float dt, int t, 
                       float* txx, float* tzz, float* txz, float* vx, float* vz,
                       float* rho, float* vs, float* vp, float* wavelet,
                       int sIdx, int sIdz);

float*
fd(int nt, int nx, int nz, float dt, float dx, float dz, 
   float* vp, float* vs, float* rho, float* wavelet, int sIdx, int sIdz) 
{
    size_t n = (size_t)nx * (size_t)nz;

    /* initialize stagerred grid arrays */
    float* txx = (float*)malloc(n * sizeof(float));
    float* tzz = (float*)malloc(n * sizeof(float));
    float* txz = (float*)malloc(n * sizeof(float));
    float* vx  = (float*)malloc(n * sizeof(float));
    float* vz  = (float*)malloc(n * sizeof(float));

    for (int t = 0; t < nt; t++) {
        fd2d_2E2T(nx, nz, dx, dz, dt, t, txx, tzz, txz, vx, vz,
                 rho, vs, vp, wavelet, sIdx, sIdz);
    }

    free(tzz);
    free(txz);
    free(vx);
    free(vz);

    return txx;
}

static void 
fd2d_2E2T(int nx, int nz, float dx, float dz, float dt, int t, 
          float* txx, float* tzz, float* txz, float* vx,
          float* vz, float* rho, float* vs, float* vp, float* wavelet,
          int sIdz, int sIdx)
{
    for (int i = 0; i < nz - 2; i++) {
        for (int j = 0; j < nx - 2; j++) {

            int idx    = i     + j     * nz;  // center
            int idx_xm = i     + (j-1) * nz;  // x-1
            int idx_xp = i     + (j+1) * nz;  // x+1
            int idx_zm = (i-1) + j     * nz;  // z-1
            int idx_zp = (i+1) + j     * nz;  // z+1

            /* applying source */
            txx[sIdz + sIdx * nz] += wavelet[t] / (dx * dx);
            tzz[sIdz + sIdx * nz] += wavelet[t] / (dx * dx);

            /* 2E2T velocity settings */
            float d_txx_dx = (txx[idx] - txx[idx_xm]) / dx;
            float d_tzz_dz = (tzz[idx] - tzz[idx_xm]) / dz; 
            float d_txz_dz = (txz[idx_zp] - txz[idx]) / dz;

            float b_x = 0.5f * (rho[idx_xp] + rho[idx]);
            float b_z = 0.5f * (rho[idx_zp] + rho[idx]);

            vx[idx] = b_x * (d_txx_dx + d_txz_dz); 
            vz[idx] = b_z * (d_txz_dz + d_tzz_dz);

            /* 2E2T pressure settings */
            float d_vx_dx = (vx[idx_xp] - vx[idx]) / dx;
            float d_vx_dz = (vx[idx_zp] - vx[idx]) / dz;
            float d_vz_dz = (vz[idx_zp] - vz[idx]) / dz;
            float d_vz_dx = (vz[idx_xp] - vz[idx]) / dx;
            
            /* elastic parameters settings */
            float lambda_xx = (vp[idx] * vp[idx] - 2.0f * vs[idx] * vs[idx]) * rho[idx];
            float lambda_zz = lambda_xx;
            float mi_xx     = (vs[idx] * vs[idx]) * rho[idx];
            float mi_zz     = mi_xx;
            float mi_vx     = 0.5f * (vs[idx]*vs[idx] * rho[idx]) + (vs[idx_xp]*vs[idx_xp] * rho[idx_xp]);
            float mi_vz     = 0.5f * (vs[idx]*vs[idx] * rho[idx]) + (vs[idx_zp]*vs[idx_zp] * rho[idx_zp]);
            float mi_xz     = powf(0.25f * (1/mi_xx + 1/mi_zz + 1/mi_vx + 1/mi_vz), -1.0f);

            txx[idx] = (lambda_xx + 2.0f * mi_xx * d_vx_dx) + (lambda_xx * d_vz_dz);
            tzz[idx] = (lambda_zz + 2.0f * mi_xx * d_vz_dz) + (lambda_xx * d_vx_dx);
            txz[idx] = mi_xz * (d_vx_dz + d_vz_dx);
        }
    }
}
