#include <stdlib.h>

#include "fd.h"
#include "bin.h"
#include "par.h"

#define SAFE_INV(x, eps)   (1.0f / ((x) + (eps)))

#define OUTPUT_PATH        "data/output/snapshots/"

#define PRINT(x)           printf("%f\n", (x))

#define IDX(i,j)           ((i) + (j) * nz)

#define AVG(a, b)          (0.5f * ((a) + (b)))

static void 
get_snapshots
( int    snap_ratio, 
  float* pressure, 
  int    time_step, 
  int    nx, 
  int    nz ) 
{ 
    if (!(time_step % snap_ratio)) {
        write_f32_bin_model(OUTPUT_PATH, pressure, nx, nz);
    }
}

static void 
fd_velocity_2E2T
( modelPar* model, 
  fdFields* fld, 
  int       nx, 
  int       nz, 
  float     dx, 
  float     dz, 
  float     dt )
{
    for (int i = 1; i < nx - 2; i++) {
        for (int j = 1; j < nz - 2; j++) {
            int idx    = IDX(i  , j  );
            int idx_xm = IDX(i  , j-1);
            int idx_xp = IDX(i  , j+1);
            int idx_zm = IDX(i  , j-1);
            int idx_zp = IDX(i-1, j  );

            float d_txx_dx = (fld->txx[idx] - fld->txx[idx_xm]) / dx;
            float d_tzz_dz = (fld->tzz[idx] - fld->tzz[idx_zm]) / dz;  
            float d_txz_dz = (fld->txz[idx_zp] - fld->txz[idx]) / dz;
            float d_txz_dx = (fld->txz[idx_xp] - fld->txz[idx]) / dx;

            float b_x = AVG(model->rho[idx_xp], model->rho[idx]);
            float b_z = AVG(model->rho[idx_zp], model->rho[idx]);
            
            fld->vx[idx] += dt * b_x *
                            (d_txx_dx + d_txz_dz);

            fld->vz[idx] += dt * b_z * 
                            (d_txz_dx + d_tzz_dz);
        }
    }
}

static void 
fd_pressure_2E2T
( modelPar* model, 
  fdFields* fld, 
  int       nx, 
  int       nz, 
  float     dx, 
  float     dz, 
  float     dt )
{
    for (int i = 1; i < nx - 2; i++)
    {
        for (int j = 1; j < nz - 2; j++)
        {
            int idx    = IDX(i  , j  );
            int idx_xp = IDX(i  , j+1);
            int idx_zp = IDX(i+1, j  );

            float d_vx_dx = (fld->vx[idx_xp] - fld->vx[idx]) / dx;
            float d_vx_dz = (fld->vx[idx_zp] - fld->vx[idx]) / dz;
            float d_vz_dz = (fld->vz[idx_zp] - fld->vz[idx]) / dz;
            float d_vz_dx = (fld->vz[idx_xp] - fld->vz[idx]) / dx;

            float lambda_xx = (model->vp[idx] * model->vp[idx] -
                               2.0f * model->vs[idx] * model->vs[idx]) *
                               model->rho[idx];

            float lambda_zz = lambda_xx;

            float mi_xx = (model->vs[idx] * model->vs[idx]) * model->rho[idx];
            float mi_zz = mi_xx;

            float mi_vx = AVG(
                model->vs[idx]   * model->vs[idx]   * model->rho[idx],
                model->vs[idx_xp] * model->vs[idx_xp] * model->rho[idx_xp]
            );

            float mi_vz = AVG(
                model->vs[idx]   * model->vs[idx]   * model->rho[idx],
                model->vs[idx_zp] * model->vs[idx_zp] * model->rho[idx_zp]
            );

            float mi_xz = SAFE_INV(0.25f *
                                   (1/mi_xx + 1/mi_zz + 1/mi_vx + 1/mi_vz),
                                   1e-8f);

            fld->txx[idx] += dt * ((lambda_xx + 2.0f * mi_xx) * d_vx_dx) +
                             (lambda_xx * d_vz_dz);

            fld->tzz[idx] += dt * ((lambda_zz + 2.0f * mi_xx) * d_vz_dz) +
                             (lambda_xx * d_vx_dx);

            fld->txz[idx] += dt * mi_xz *
                             (d_vx_dz + d_vz_dx);
        }
    }
}

float* 
fd
( fdFields*   fld, 
  modelPar*   model, 
  geomPar*    geom, 
  snapshots*  snap, 
  waveletPar* wav )
{
    size_t n = model->nx * model->nz;

    /* initialize staggered grid arrays */
    fld->txx = (float*)calloc(n, sizeof(float));
    fld->tzz = (float*)calloc(n, sizeof(float));
    fld->txz = (float*)calloc(n, sizeof(float));
    fld->vx  = (float*)calloc(n, sizeof(float));
    fld->vz  = (float*)calloc(n, sizeof(float));

    if (!fld->txx || !fld->tzz || !fld->txz || !fld->vx || !fld->vz) {
        return NULL;
    }

    int nx = model->nx;
    int nz = model->nz;
    float dx = model->dx;
    float dz = model->dz;

    int snap_ratio = wav->nt / snap->snap_num;

    for (int t = 0; t < wav->nt; t++) {
        int s_idx = geom->sIdz + geom->sIdx*nz;

        fld->txx[s_idx] += wav->wavelet[t] / (dx*dz);
        fld->tzz[s_idx] += wav->wavelet[t] / (dx*dz);

        fd_velocity_2E2T(model, fld, nx, nz, dx, dz, wav->dt);
        fd_pressure_2E2T(model, fld, nx, nz, dx, dz, wav->dt);

        if (snap->snap_bool) {
            get_snapshots(snap_ratio, fld->txx, t, nx, nz);
        }
    }

    free(fld->tzz);
    free(fld->txz);
    free(fld->vx);
    free(fld->vz);

    return fld->txx; 
}

