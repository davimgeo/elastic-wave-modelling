#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "io_bin.h"
#include "model_ext.h"
#include "fd_damp.h"
#include "fd_core.h"
#include "cfg_par.h"

/* finite-difference coefficients */
#define FDM8E1 6.97545e-4f
#define FDM8E2 9.57031e-3f
#define FDM8E3 7.97526e-2f
#define FDM8E4 1.19628906f

#define OUTPUT_PATH "data/output/snapshots/"

static void generate_p(const fields_t *fld, int nxx, int nzz)
{
  for (int j = 0; j < nxx; j++)
  {
    for (int i = 0; i < nzz; i++)
    {
      fld->calc_p[i + j * nzz] =
        0.5f * (fld->txx[i + j * nzz] + fld->tzz[i + j * nzz]);
    }
  }
}

static void get_snapshots(const config_t *cfg, const fields_t *fld, int time_step)
{
  if (!(time_step % cfg->snap_ratio))
  {
  generate_p(fld, cfg->nxx, cfg->nzz);

  const char *filenames[] = {
    "p_%dx%d_tid_%d.bin",
    "vx_%dx%d_tid_%d.bin",
    "vz_%dx%d_tid_%d.bin"
  };

  float *fields[] = { fld->calc_p, fld->vx, fld->vz };

  char current_snap[256], full_path[512];

  #define FIELDS_SIZE 3
  for (int i = 0; i < FIELDS_SIZE; i++)
  {
    snprintf(current_snap, sizeof(current_snap), filenames[i],
              cfg->nxx, cfg->nzz, time_step);
    snprintf(full_path, sizeof(full_path), "%s%s", OUTPUT_PATH, current_snap);

    write2D(full_path, fields[i], sizeof(float), cfg->nzz, cfg->nxx);
  }

  printf("Generating snapshot num %d...\n", time_step);
  }
}

static inline void
fd_velocity8E2T
( fields_t *restrict    fld,
  const model_t *restrict m,
  int                   nxx,
  int                   nzz,
  float              inv_dx,
  float              inv_dz,
  float                  dt,
  damping_t           *damp )
{
  float *restrict txx = fld->txx;
  float *restrict tzz = fld->tzz;
  float *restrict txz = fld->txz;
  float *restrict vx  = fld->vx;
  float *restrict vz  = fld->vz;

  const float *restrict rho = m->rho;

  #pragma omp for schedule(static) nowait
  for (int j = 4; j < nxx - 4; j++)
  {
    for (int i = 4; i < nzz - 4; i++)
    {
      float dtxx_dx =
        (FDM8E1 * (txx[i + (j - 4) * nzz] - txx[i + (j + 3) * nzz]) +
         FDM8E2 * (txx[i + (j + 2) * nzz] - txx[i + (j - 3) * nzz]) +
         FDM8E3 * (txx[i + (j - 2) * nzz] - txx[i + (j + 1) * nzz]) +
         FDM8E4 * (txx[i + j * nzz] - txx[i + (j - 1) * nzz])) * inv_dx;

      float dtxz_dz =
        (FDM8E1 * (txz[(i - 3) + j * nzz] - txz[(i + 4) + j * nzz]) +
         FDM8E2 * (txz[(i + 3) + j * nzz] - txz[(i - 2) + j * nzz]) +
         FDM8E3 * (txz[(i - 1) + j * nzz] - txz[(i + 2) + j * nzz]) +
         FDM8E4 * (txz[(i + 1) + j * nzz] - txz[i + j * nzz])) * inv_dz;

      float dtxz_dx =
        (FDM8E1 * (txz[i + (j - 3) * nzz] - txz[i + (j + 4) * nzz]) +
         FDM8E2 * (txz[i + (j + 3) * nzz] - txz[i + (j - 2) * nzz]) +
         FDM8E3 * (txz[i + (j - 1) * nzz] - txz[i + (j + 2) * nzz]) +
         FDM8E4 * (txz[i + (j + 1) * nzz] - txz[i + j * nzz])) * inv_dx;

      float dtzz_dz =
        (FDM8E1 * (tzz[(i - 4) + j * nzz] - tzz[(i + 3) + j * nzz]) +
         FDM8E2 * (tzz[(i + 2) + j * nzz] - tzz[(i - 3) + j * nzz]) +
         FDM8E3 * (tzz[(i - 2) + j * nzz] - tzz[(i + 1) + j * nzz]) +
         FDM8E4 * (tzz[i + j * nzz] - tzz[(i - 1) + j * nzz])) * inv_dz;

      float rho_inv  = 1.0f / (0.5f * (rho[i + j * nzz] + rho[i + (j + 1) * nzz]));
      float rho_inv2 = 1.0f / (0.5f * (rho[i + j * nzz] + rho[(i + 1) + j * nzz]));

      vx[i + j * nzz] += dt * rho_inv  * (dtxx_dx + dtxz_dz);
      vz[i + j * nzz] += dt * rho_inv2 * (dtxz_dx + dtzz_dz);

      vx[i + j * nzz] *= (damp->x[j] * damp->z[i]);
      vz[i + j * nzz] *= (damp->x[j] * damp->z[i]);
    }
  }
}

static inline void
fd_pressure8E2T
( fields_t *restrict    fld,
  const model_t *restrict m,
  int                   nxx,
  int                   nzz,
  float              inv_dx,
  float              inv_dz,
  float                  dt,
  damping_t           *damp )
{
  float *restrict txx = fld->txx;
  float *restrict tzz = fld->tzz;
  float *restrict txz = fld->txz;
  float *restrict vx  = fld->vx;
  float *restrict vz  = fld->vz;

  const float *restrict vp  = m->vp;
  const float *restrict vs  = m->vs;
  const float *restrict rho = m->rho;

  #pragma omp for schedule(static)
  for (int j = 4; j < nxx - 4; j++)
  {
    for (int i = 4; i < nzz - 4; i++)
    {
      float dvx_dx =
        (FDM8E1 * (vx[i + (j - 3) * nzz] - vx[i + (j + 4) * nzz]) +
         FDM8E2 * (vx[i + (j + 3) * nzz] - vx[i + (j - 2) * nzz]) +
         FDM8E3 * (vx[i + (j - 1) * nzz] - vx[i + (j + 2) * nzz]) +
         FDM8E4 * (vx[i + (j + 1) * nzz] - vx[i + j * nzz])) * inv_dx;

      float dvz_dz =
        (FDM8E1 * (vz[(i - 3) + j * nzz] - vz[(i + 4) + j * nzz]) +
         FDM8E2 * (vz[(i + 3) + j * nzz] - vz[(i - 2) + j * nzz]) +
         FDM8E3 * (vz[(i - 1) + j * nzz] - vz[(i + 2) + j * nzz]) +
         FDM8E4 * (vz[(i + 1) + j * nzz] - vz[i + j * nzz])) * inv_dz;

      float dvx_dz =
        (FDM8E1 * (vx[(i - 4) + j * nzz] - vx[(i + 3) + j * nzz]) +
         FDM8E2 * (vx[(i + 2) + j * nzz] - vx[(i - 3) + j * nzz]) +
         FDM8E3 * (vx[(i - 2) + j * nzz] - vx[(i + 1) + j * nzz]) +
         FDM8E4 * (vx[i + j * nzz] - vx[(i - 1) + j * nzz])) * inv_dz;

      float dvz_dx =
        (FDM8E1 * (vz[i + (j - 4) * nzz] - vz[i + (j + 3) * nzz]) +
         FDM8E2 * (vz[i + (j + 2) * nzz] - vz[i + (j - 3) * nzz]) +
         FDM8E3 * (vz[i + (j - 2) * nzz] - vz[i + (j + 1) * nzz]) +
         FDM8E4 * (vz[i + j * nzz] - vz[i + (j - 1) * nzz])) * inv_dx;

      float vp2       = vp[i + j * nzz] * vp[i + j * nzz];
      float vs2       = vs[i + j * nzz] * vs[i + j * nzz];
      float vs2_xp    = vs[(i + 1) + j * nzz] * vs[(i + 1) + j * nzz];
      float vs2_zp    = vs[i + (j + 1) * nzz] * vs[i + (j + 1) * nzz];
      float vs2_xp_zp = vs[(i + 1) + (j + 1) * nzz] * vs[(i + 1) + (j + 1) * nzz];

      float lamb = rho[i + j * nzz] * (vp2 - 2.0f * vs2);
      float mi   = rho[i + j * nzz] * vs2;

      float mi1    = rho[i + j * nzz]         * vs2;
      float mi2    = rho[(i + 1) + j * nzz]   * vs2_xp;
      float mi3    = rho[i + (j + 1) * nzz]   * vs2_zp;
      float mi4    = rho[(i + 1) + (j + 1) * nzz] * vs2_xp_zp;
      float mi_avg = 4.0f / ((1.0f / mi1) + (1.0f / mi2) + (1.0f / mi3) + (1.0f / mi4));

      txx[i + j * nzz] += dt * ((lamb + 2.0f * mi) * dvx_dx + lamb * dvz_dz);
      tzz[i + j * nzz] += dt * ((lamb + 2.0f * mi) * dvz_dz + lamb * dvx_dx);
      txz[i + j * nzz] += dt * mi_avg * (dvx_dz + dvz_dx);

      txx[i + j * nzz] *= (damp->x[j] * damp->z[i]);
      tzz[i + j * nzz] *= (damp->x[j] * damp->z[i]);
      txz[i + j * nzz] *= (damp->x[j] * damp->z[i]);
    }
  }
}

static void inject_source
( float              sIdx_f,
  float              sIdz_f,
  float              inv_dx_dz,
  const config_t    *cfg,
  fields_t          *fld,
  size_t             t )
{
  int sIdx = (int)sIdx_f;
  int sIdz = (int)sIdz_f;

  int s_idx = (sIdz + cfg->nb) + (sIdx + cfg->nb) * cfg->nzz;

  fld->txx[s_idx] += cfg->wavelet[t] * inv_dx_dz;
  fld->tzz[s_idx] += cfg->wavelet[t] * inv_dx_dz;
}

static void init_fields(fields_t *fld, int nxx, int nzz)
{
  size_t n = (size_t)nxx * (size_t)nzz;

  fld->txx = (float *)calloc(n, sizeof(float));
  fld->tzz = (float *)calloc(n, sizeof(float));
  fld->txz = (float *)calloc(n, sizeof(float));
  fld->vx  = (float *)calloc(n, sizeof(float));
  fld->vz  = (float *)calloc(n, sizeof(float));

  fld->calc_p = (float *)calloc(n, sizeof(float));

  if (!fld->txx || !fld->tzz || !fld->txz || !fld->vx || !fld->vz || !fld->calc_p) 
  {
    perror("Could not allocate fields\n");
    exit(EXIT_FAILURE);
  }
}

static void free_fields(fields_t *fld)
{
  free(fld->txx);
  free(fld->tzz);
  free(fld->txz);
  free(fld->vx);
  free(fld->vz);
  free(fld->calc_p);
}

static void
register_seismogram
( int t,
  const config_t *cfg,
  float  inv_dx_dz,
  float *seismogram,
  float *field )
{
  #pragma omp for schedule(static)
  for (int i = 0; i < cfg->r_f_lines - 1; i++)
  {
    int ix = (int)(cfg->rcv_x[i] * inv_dx_dz) + cfg->nb;
    int iz = (int)(cfg->rcv_z[i] * inv_dx_dz) + cfg->nb;

    seismogram[t + i * cfg->nt] = field[iz + ix * cfg->nzz];
  }
}

void fd(const config_t *cfg, fields_t *fld, model_t *m)
{
  set_boundary(cfg, m);

  init_fields(fld, cfg->nxx, cfg->nzz);

  int seismogram_size = cfg->nt * cfg->r_f_lines;
  float* seismogram = (float *)calloc((size_t)seismogram_size, sizeof(float));
  if (!seismogram) 
  {
    perror("Could not allocate Seismogram\n");
    exit(EXIT_FAILURE);
  }

  damping_t damp = get_damp(cfg);

  float inv_dx = 1.0f / cfg->dx;
  float inv_dz = 1.0f / cfg->dz;
  float inv_dx_dz = inv_dx * inv_dz;

  for (size_t i = 0; i < cfg->src_f_lines - 1; i++)
  {
    float sIdx = cfg->src_x[i];
    float sIdz = cfg->src_z[i];

    #pragma omp parallel
    {
      for (size_t t = 0; t < (size_t)cfg->nt; t++)
      {
        #pragma omp single
        inject_source(sIdx, sIdz, inv_dx_dz, cfg, fld, t);

        fd_pressure8E2T(fld, m, cfg->nxx, cfg->nzz, inv_dx, inv_dz, cfg->dt, &damp);
        fd_velocity8E2T(fld, m, cfg->nxx, cfg->nzz, inv_dx, inv_dz, cfg->dt, &damp);

        register_seismogram(t, cfg, inv_dx_dz, seismogram, fld->vx);

        #pragma omp single
        if (cfg->snap_bool)
          get_snapshots(cfg, fld, (int)t);
      }
    }
  }

  write2D("data/output/seismogram_vx_1150x648.bin", seismogram, sizeof(float), cfg->nt, cfg->r_f_lines);

  free(seismogram);
  free(damp.x);
  free(damp.z);
  free_fields(fld);
}

