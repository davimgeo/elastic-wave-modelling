#include <stdlib.h>

#include "fd.h"
#include "par.h"

static void 
get_snapshots(snapshots* snap, float* pressure, int time_step)
{
    if (snap->snap_bool && !(time_step % snap->snap_ratio)) {
        // save pressure on data/snapshots
        // using write_f32_bin_model() function
    }
}

static void  
fd2d_2E2T(modelPar* model, geomPar* geom, waveletPar* wav, fdFields* fld, int time_step)
{
    int nz = model->nz;
    int nx = model->nx;
    float dx = model->dx;
    float dz = model->dz;

    for (int i = 0; i < model->nz - 2; i++) {
        for (int j = 0; j < model->nx - 2; j++) {
            int idx    = i     + j     * nz;  // center
            int idx_xm = i     + (j-1) * nz;  // x-1
            int idx_xp = i     + (j+1) * nz;  // x+1
            int idx_zm = (i-1) + j     * nz;  // z-1
            int idx_zp = (i+1) + j     * nz;  // z+1

            /* applying source */
            fld->txx[geom->sIdz + geom->sIdx*nz] += wav->wavelet[time_step] / (dx * dx);
            fld->tzz[geom->sIdz + geom->sIdx*nz] += wav->wavelet[time_step] / (dx * dx);

            /* 2E2T velocity settings */
            float d_txx_dx = (fld->txx[idx] - fld->txx[idx_xm]) / (2.0f*dx);
            float d_tzz_dz = (fld->tzz[idx] - fld->tzz[idx_xm]) / (2.0f*dz);
            float d_txz_dz = (fld->txz[idx_zp] - fld->txz[idx]) / (2.0f*dz);

            float b_x = 0.5f * (model->rho[idx_xp] + model->rho[idx]);
            float b_z = 0.5f * (model->rho[idx_zp] + model->rho[idx]);

            fld->vx[idx] = b_x * (d_txx_dx + d_txz_dz); 
            fld->vz[idx] = b_z * (d_txz_dz + d_tzz_dz);

            /* 2E2T pressure settings */
            float d_vx_dx = (fld->vx[idx_xp] - fld->vx[idx]) / dx;
            float d_vx_dz = (fld->vx[idx_zp] - fld->vx[idx]) / dz;
            float d_vz_dz = (fld->vz[idx_zp] - fld->vz[idx]) / dz;
            float d_vz_dx = (fld->vz[idx_xp] - fld->vz[idx]) / dx;
            
            /* elastic parameters settings */
            float lambda_xx = (model->vp[idx] * model->vp[idx] - 2.0f * model->vs[idx] * model->vs[idx]) * model->rho[idx];
            float lambda_zz = lambda_xx;
            float mi_xx     = (model->vs[idx] * model->vs[idx]) * model->rho[idx];
            float mi_zz     = mi_xx;
            float mi_vx     = 0.5f * ((model->vs[idx]    * model->vs[idx]    * model->rho[idx]) +
                                      (model->vs[idx_xp] * model->vs[idx_xp] * model->rho[idx_xp]));
            float mi_vz     = 0.5f * ((model->vs[idx]    * model->vs[idx]    * model->rho[idx]) +
                                      (model->vs[idx_zp] * model->vs[idx_zp] * model->rho[idx_zp]));
            float mi_xz = 1.0f / (0.25f * (1.0f/mi_xx + 1.0f/mi_zz + 1.0f/mi_vx + 1.0f/mi_vz));

            fld->txx[idx] = (lambda_xx + 2.0f * mi_xx * d_vx_dx) + (lambda_xx * d_vz_dz);
            fld->tzz[idx] = (lambda_zz + 2.0f * mi_xx * d_vz_dz) + (lambda_xx * d_vx_dx);
            fld->txz[idx] = mi_xz * (d_vx_dz + d_vz_dx);
        }
    }
}

float* 
fd(modelPar* model, geomPar* geom, snapshots* snap, waveletPar* wav)
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
