#include <stdlib.h>
#include <math.h>

#include "fd.h"
#include "par.h"

static void get_snapshots(snapshots* snap, float* pressure, int time_step)
{
    if (snap->snap_bool && !(time_step % snap->snap_ratio)) {
        // save pressure on data/snapshots
    }
}

static void  
fd2d_2E2T(modelPar* model, geomPar* geom, waveletPar* wav, fdFields* fld, int time_step)
{
    int nx = model->nx;
    int nz = model->nz;
    float dx = model->dx;
    float dz = model->dz;
    float dt = wav->dt;

    float *txx = fld->txx;
    float *tzz = fld->tzz;
    float *txz = fld->txz;
    float *vx  = fld->vx;
    float *vz  = fld->vz;

    float *rho = model->rho;
    float *vs  = model->vs;
    float *vp  = model->vp;

    float *wavelet = wav->wavelet;

    int sIdx = geom->sIdx;
    int sIdz = geom->sIdz;

    for (int i = 0; i < nz - 2; i++) {
        for (int j = 0; j < nx - 2; j++) {

            int idx    = i     + j     * nz;  // center
            int idx_xm = i     + (j-1) * nz;  // x-1
            int idx_xp = i     + (j+1) * nz;  // x+1
            int idx_zm = (i-1) + j     * nz;  // z-1
            int idx_zp = (i+1) + j     * nz;  // z+1

            /* applying source */
            txx[sIdz + sIdx*nz] += wavelet[time_step] / (dx * dx);
            tzz[sIdz + sIdx*nz] += wavelet[time_step] / (dx * dx);

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
            
            // elastic parameters settings
            float lambda_xx = (vp[idx] * vp[idx] - 2.0f * vs[idx] * vs[idx]) * rho[idx];
            float lambda_zz = lambda_xx;
            float mi_xx     = (vs[idx] * vs[idx]) * rho[idx];
            float mi_zz     = mi_xx;
            float mi_vx     = 0.5f * (vs[idx]*vs[idx] * rho[idx]) +
                              (vs[idx_xp]*vs[idx_xp] * rho[idx_xp]);
            float mi_vz     = 0.5f * (vs[idx]*vs[idx] * rho[idx]) +
                              (vs[idx_zp]*vs[idx_zp] * rho[idx_zp]);
            float mi_xz     = powf(0.25f * (1/mi_xx + 1/mi_zz + 1/mi_vx + 1/mi_vz), -1.0f);

            txx[idx] = (lambda_xx + 2.0f * mi_xx * d_vx_dx) + (lambda_xx * d_vz_dz);
            tzz[idx] = (lambda_zz + 2.0f * mi_xx * d_vz_dz) + (lambda_xx * d_vx_dx);
            txz[idx] = mi_xz * (d_vx_dz + d_vz_dx);
        }
    }
}

float* fd(modelPar* model, geomPar* geom, snapshots* snap, waveletPar* wav)
{
    size_t n = (size_t)model->nx * (size_t)model->nz;

    fdFields fld;
    /* initialize staggered grid arrays */
    fld.txx = (float*)calloc(n, sizeof(float));
    fld.tzz = (float*)calloc(n, sizeof(float));
    fld.txz = (float*)calloc(n, sizeof(float));
    fld.vx  = (float*)calloc(n, sizeof(float));
    fld.vz  = (float*)calloc(n, sizeof(float));

    if (fld.txx == NULL || fld.tzz == NULL || fld.txz == NULL ||
        fld.vx  == NULL || fld.vz  == NULL) {
        return NULL;
    }

    for (int t = 0; t < wav->nt; t++) {
        fd2d_2E2T(model, geom, wav, &fld, t);

        get_snapshots(snap, fld.txx, t);
    }

    free(fld.tzz);
    free(fld.txz);
    free(fld.vx);
    free(fld.vz);

    return fld.txx; /* Caller is responsible for freeing txx */
}
