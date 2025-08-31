#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "bin.h"
#include "fd.h"
#include "par.h"

#define PI 3.1415926f

/* Finite-difference coefficients */
#define FDM8E1 6.97545e-4f
#define FDM8E2 9.57031e-3f
#define FDM8E3 7.97526e-2f
#define FDM8E4 1.19628906f

#define OUTPUT_PATH "data/output/snapshots/"

#define PRINT(x) printf("%f\n", (x))

#define AVG(a, b) (0.5f * ((a) + (b)))

#define DOUBLE(x) ((x) * (x))

float* ricker(int nt, float dt, float fmax) 
{
  float* ricker = (float *)calloc((size_t)nt, sizeof(float));
  if (ricker == NULL) 
    return NULL;
  
  float t0 = 2.0f * PI / fmax;
  float fc = fmax / (3.0f*sqrtf(PI));

  for (int i = 0; i < nt; i++) 
  {
    float t = (i * dt) - t0;
    float arg = PI * (PI * PI * fc * fc * t * t);
    ricker[i] = (1.0f - 2.0f * arg) * expf(-arg);
  }
  return ricker;
}

static void 
get_snapshots
( int    snap_ratio, 
  float *pressure, 
  int    time_step,
  int    nx, 
  int    nz )
{
  char current_snap[256];
  snprintf(current_snap, sizeof(current_snap), 
      "txx_%dx%d_tid_%d.bin", nx, nz, time_step);

  if (!(time_step % snap_ratio))
  {
    char full_path[512];
    snprintf(full_path, sizeof(full_path), "%s%s", OUTPUT_PATH, current_snap);

    write_f32_bin_model(full_path, pressure, nx, nz);
    printf("Generating snapshot num %d...\n", time_step);
  }
}

static void 
fd_velocity_8E2T
( float    *rho, 
  fdFields *fld, 
  int       nx, 
  int       nz,
  float     dx, 
  float     dz, 
  float     dt )
{
  #pragma omp parallel for
  for (int index = 0; index < nx * nz; index++)
  {
    int i = (int)(index % nz);
    int j = (int)(index / nz);

    if ((i >= 4) && (i < nz - 4) && (j > 4) && (j < nx - 4))
    {
      float dtxx_dx =
          (FDM8E1 * (fld->txx[i + (j - 4) * nz] - fld->txx[i + (j + 3) * nz]) +
           FDM8E2 * (fld->txx[i + (j + 2) * nz] - fld->txx[i + (j - 3) * nz]) +
           FDM8E3 * (fld->txx[i + (j - 2) * nz] - fld->txx[i + (j + 1) * nz]) +
           FDM8E4 * (fld->txx[i + j * nz] - fld->txx[i + (j - 1) * nz])) / dx;

      float dtxz_dz =
          (FDM8E1 * (fld->txz[(i - 3) + j * nz] - fld->txz[(i + 4) + j * nz]) +
           FDM8E2 * (fld->txz[(i + 3) + j * nz] - fld->txz[(i - 2) + j * nz]) +
           FDM8E3 * (fld->txz[(i - 1) + j * nz] - fld->txz[(i + 2) + j * nz]) +
           FDM8E4 * (fld->txz[(i + 1) + j * nz] - fld->txz[i + j * nz])) / dz;

      float inv_rho_1 = 1.0f / AVG(rho[i + j * nz], rho[i + (j+1) * nz]);

      fld->vx[i + j * nz] += dt * inv_rho_1 * (dtxx_dx + dtxz_dz);
    }

    if ((i > 4) && (i < nz - 4) && (j >= 4) && (j < nx - 4))
    {
      float dtxz_dx =
          (FDM8E1 * (fld->txz[i + (j - 3) * nz] - fld->txz[i + (j + 4) * nz]) +
           FDM8E2 * (fld->txz[i + (j + 3) * nz] - fld->txz[i + (j - 2) * nz]) +
           FDM8E3 * (fld->txz[i + (j - 1) * nz] - fld->txz[i + (j + 2) * nz]) +
           FDM8E4 * (fld->txz[i + (j + 1) * nz] - fld->txz[i + j * nz])) / dx;

      float dtzz_dz =
          (FDM8E1 * (fld->tzz[(i - 4) + j * nz] - fld->tzz[(i + 3) + j * nz]) +
           FDM8E2 * (fld->tzz[(i + 2) + j * nz] - fld->tzz[(i - 3) + j * nz]) +
           FDM8E3 * (fld->tzz[(i - 2) + j * nz] - fld->tzz[(i + 1) + j * nz]) +
           FDM8E4 * (fld->tzz[i + j * nz] - fld->tzz[(i - 1) + j * nz])) / dz;

      float inv_rho_2 = 1.0f / AVG(rho[i + j * nz], rho[(i+1) + j * nz]);

      fld->vz[i + j * nz] += dt * inv_rho_2 * (dtxz_dx + dtzz_dz);
    }
  }
}

static void 
fd_pressure_8E2T
( float    *vp,
  float    *vs,
  float    *rho,
  fdFields *fld, 
  int       nx, 
  int       nz,
  float     dx, 
  float     dz, 
  float     dt )
{
  #pragma omp parallel for
  for (int index = 0; index < nx * nz; index++)
  {
    int i = (int)(index % nz);
    int j = (int)(index / nz);

    if ((i >= 4) && (i < nz - 4) && (j >= 4) && (j < nx - 4))
    {
      float dvx_dx =
          (FDM8E1 * (fld->vx[i + (j - 3) * nz] - fld->vx[i + (j + 4) * nz]) +
           FDM8E2 * (fld->vx[i + (j + 3) * nz] - fld->vx[i + (j - 2) * nz]) +
           FDM8E3 * (fld->vx[i + (j - 1) * nz] - fld->vx[i + (j + 2) * nz]) +
           FDM8E4 * (fld->vx[i + (j + 1) * nz] - fld->vx[i + j * nz])) / dx;

      float dvz_dz =
          (FDM8E1 * (fld->vz[(i - 3) + j * nz] - fld->vz[(i + 4) + j * nz]) +
           FDM8E2 * (fld->vz[(i + 3) + j * nz] - fld->vz[(i - 2) + j * nz]) +
           FDM8E3 * (fld->vz[(i - 1) + j * nz] - fld->vz[(i + 2) + j * nz]) +
           FDM8E4 * (fld->vz[(i + 1) + j * nz] - fld->vz[i + j * nz])) / dz;

      float vp2 = vp[i + j * nz] * vp[i + j * nz];
      float vs2 = vs[i + j * nz] * vs[i + j * nz];

      float lambda = rho[i + j * nz] * (vp2 - 2.0f * vs2);

      float mi     = rho[i + j * nz] * vs2;

      fld->txx[i + j * nz] +=
          dt * ((lambda + 2.0f * mi) * dvx_dx +
                 lambda * dvz_dz);

      fld->tzz[i + j * nz] +=
          dt * ((lambda + 2.0f * mi) * dvz_dz +
                 lambda * dvx_dx);
    }

    if ((i >= 4) && (i <= nz - 4) && (j >= 4) && (j <= nx - 4))
    {
      float dvx_dz =
          (FDM8E1 * (fld->vx[(i - 4) + j * nz] - fld->vx[(i + 3) + j * nz]) +
           FDM8E2 * (fld->vx[(i + 2) + j * nz] - fld->vx[(i - 3) + j * nz]) +
           FDM8E3 * (fld->vx[(i - 2) + j * nz] - fld->vx[(i + 1) + j * nz]) +
           FDM8E4 * (fld->vx[i + j * nz] - fld->vx[(i - 1) + j * nz])) / dz;

      float dvz_dx =
          (FDM8E1 * (fld->vz[i + (j - 4) * nz] - fld->vz[i + (j + 3) * nz]) +
           FDM8E2 * (fld->vz[i + (j + 2) * nz] - fld->vz[i + (j - 3) * nz]) +
           FDM8E3 * (fld->vz[i + (j - 2) * nz] - fld->vz[i + (j + 1) * nz]) +
           FDM8E4 * (fld->vz[i + j * nz] - fld->vz[i + (j - 1) * nz])) / dx;

      float vs2       = vs[i + j * nz] * vs[i + j * nz];
      float vs2_xp    = vs[(i+1) + j * nz] * vs[i + j * nz];
      float vs2_zp    = vs[i + (j+1) * nz] * vs[i + (j+1) * nz];
      float vs2_xp_zp = vs[(i+1) + (j+1) * nz] * vs[(i+1) + (j+1) * nz];

      float mi1 = rho[i + j * nz] * vs2;
      float mi2 = rho[(i+1) + j * nz] * vs2_xp;
      float mi3 = rho[i + (j+1) * nz] * vs2_zp;
      float mi4 = rho[(i+1) + (j+1) * nz] * vs2_xp_zp;

      float mi_avg = 4.0f / ((1.0f / mi1) + (1.0f / mi2) +
                             (1.0f / mi3) + (1.0f / mi4));

      fld->txz[i + j * nz] += dt * mi_avg * (dvx_dz + dvz_dx);
    }
  }
}

void 
fd
( fdFields  *fld, 
  float     *vp,
  float     *vs,
  float     *rho,
  int        nx,
  int        nz,
  int        nt,
  float     *wavelet,
  float      dt,
  float      dx,
  float      dz,
  int        sIdx,
  int        sIdz,
  snapshots *snap )
{
  size_t n = nx * nz;

  fld->txx = (float *)calloc(n, sizeof(float));
  fld->tzz = (float *)calloc(n, sizeof(float));
  fld->txz = (float *)calloc(n, sizeof(float));
  fld->vx  = (float *)calloc(n, sizeof(float));
  fld->vz  = (float *)calloc(n, sizeof(float));

  if (!fld->txx || !fld->tzz || !fld->txz || !fld->vx || !fld->vz)
  {
    perror("Could not allocate fields\n");
    return;
  }

  int snap_ratio = nt / snap->snap_num;

  for (int t = 0; t < nt; t++)
  {
    int s_idx = sIdz + sIdx * nz;

    fld->txx[s_idx] += wavelet[t] / (dx * dz);
    fld->tzz[s_idx] += wavelet[t] / (dx * dz);

    fd_velocity_8E2T(rho, fld, nx, nz, dx, dz, dt);
    fd_pressure_8E2T(vp, vs, rho, fld, nx, nz, dx, dz, dt);

    if (snap->snap_bool)
      get_snapshots(snap_ratio, fld->txx, t, nx, nz);
  }
}
