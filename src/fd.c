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

static const float vp  = 1500.0f;
static const float vs  = 0.0f;
static const float rho = 1000.0f;

static const float M = vs * vs * rho;
static const float L = vp * vp * rho - 2.0f * M;
static const float B = 1.0f / rho;

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
( modelPar *model, 
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

      fld->vx[i + j * nz] += dt * B * (dtxx_dx + dtxz_dz);
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

      fld->vz[i + j * nz] += dt * B * (dtxz_dx + dtzz_dz);
    }
  }
}

static void 
fd_pressure_8E2T
( modelPar *model, 
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

      fld->txx[i + j * nz] += dt * ((L + 2 * M) * dvx_dx + L * dvz_dz);
      fld->tzz[i + j * nz] += dt * ((L + 2 * M) * dvz_dz + L * dvx_dx);
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

      fld->txz[i + j * nz] += dt * M * (dvx_dz + dvz_dx);
    }
  }
}

void 
fd
( fdFields   *fld, 
  modelPar *model, 
  geomPar   *geom, 
  snapshots *snap,
  waveletPar *wav )
{
  size_t n = model->nx * model->nz;

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

  int nx   = model->nx;
  int nz   = model->nz;
  float dx = model->dx;
  float dz = model->dz;

  int snap_ratio = wav->nt / snap->snap_num;

  for (int t = 0; t < wav->nt; t++)
  {
    int s_idx = geom->sIdz + geom->sIdx * nz;

    fld->txx[s_idx] += wav->wavelet[t] / (dx * dz);
    fld->tzz[s_idx] += wav->wavelet[t] / (dx * dz);

    fd_velocity_8E2T(model, fld, nx, nz, dx, dz, wav->dt);
    fd_pressure_8E2T(model, fld, nx, nz, dx, dz, wav->dt);

    if (snap->snap_bool)
      get_snapshots(snap_ratio, fld->txx, t, nx, nz);
  }
}
