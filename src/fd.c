#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "bin.h"
#include "fd.h"
#include "par.h"

/* finite-difference coefficients */
#define FDM8E1 6.97545e-4f
#define FDM8E2 9.57031e-3f
#define FDM8E3 7.97526e-2f
#define FDM8E4 1.19628906f

#define OUTPUT_PATH "data/output/snapshots/"

static void set_boundary(const config_t *p, model_t *m)
{
  int nb  = p->nb;
  int nx  = p->nx;
  int nz  = p->nz;
  int nxx = p->nxx;   
  int nzz = p->nzz;  

  size_t n = nxx * nzz;

  float *vp_ext  = (float *)calloc(n, sizeof(float));
  float *vs_ext  = (float *)calloc(n, sizeof(float));
  float *rho_ext = (float *)calloc(n, sizeof(float));

  /* copy original arr into ext */
  for (int j = 0; j < nx; j++) 
  {
    for (int i = 0; i < nz; i++) 
    {
      vp_ext[(i + nb) + (j + nb) * nzz]  = m->vp[i + j * nz];
      vs_ext[(i + nb) + (j + nb) * nzz]  = m->vs[i + j * nz];
      rho_ext[(i + nb) + (j + nb) * nzz] = m->rho[i + j * nz];
    }
  }

  /* pad bottom */
  for (int j = nb; j < nx+nb; j++) 
  {
    for (int i = 0; i < nb; i++) 
    {
      vp_ext[i + j * nzz]  = vp_ext[nb + j * nzz];
      vs_ext[i + j * nzz]  = vs_ext[nb + j * nzz];
      rho_ext[i + j * nzz] = rho_ext[nb + j * nzz];

      vp_ext[(nz + nb + i) + j * nzz]  = vp_ext[(nz + nb - 1) + j * nzz];
      vs_ext[(nz + nb + i) + j * nzz]  = vs_ext[(nz + nb - 1) + j * nzz];
      rho_ext[(nz + nb + i) + j * nzz] = rho_ext[(nz + nb - 1) + j * nzz];
    }
  }

  /* pad left and right respectively */
  for (int i = 0; i < nzz; i++) 
  {
    for (int j = 0; j < nb; j++) 
    {
      // counld vectorize because of strided loop
      vp_ext[i + j * nzz]  = vp_ext[i + nb * nzz];
      vs_ext[i + j * nzz]  = vs_ext[i + nb * nzz];
      rho_ext[i + j * nzz] = rho_ext[i + nb * nzz];

      vp_ext[i + (nx + nb + j) * nzz]  = vp_ext[i + (nx + nb - 1) * nzz];
      vs_ext[i + (nx + nb + j) * nzz]  = vs_ext[i + (nx + nb - 1) * nzz];
      rho_ext[i + (nx + nb + j) * nzz] = rho_ext[i + (nx + nb - 1) * nzz];  
    }
  }

  /* swap pointers to new arr */
  free(m->vp);  m->vp  = vp_ext;
  free(m->vs);  m->vs  = vs_ext;
  free(m->rho); m->rho = rho_ext;
}

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

void get_snapshots(const config_t *p, const fields_t *fld, int time_step)
{
  if (!(time_step % p->snap_ratio)) 
  {
    generate_p(fld, p->nxx, p->nzz);

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
               p->nxx, p->nzz, time_step);
      snprintf(full_path, sizeof(full_path), "%s%s", OUTPUT_PATH, current_snap);
  
      write2D(full_path, fields[i], sizeof(float), p->nzz, p->nxx);
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

      // apply boundary
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

      // apply boundary
      txx[i + j * nzz] *= (damp->x[j] * damp->z[i]);
      tzz[i + j * nzz] *= (damp->x[j] * damp->z[i]);
      txz[i + j * nzz] *= (damp->x[j] * damp->z[i]);
    }
  }
}

void static inject_source(const config_t *p, fields_t *fld, size_t t)
{
  int s_idx = (p->sIdz + p->nb) + (p->sIdx + p->nb) * p->nzz;

  fld->txx[s_idx] += p->wavelet[t] / (p->dx * p->dz);
  fld->tzz[s_idx] += p->wavelet[t] / (p->dx * p->dz);
}

static void get_damp_x(damping_t *damp, const config_t *p)
{
  for (int j = 0; j < p->nxx; j++) 
  {
    // same as nb <= j <= nz + nb
    if ((unsigned)(j - p->nb) < (p->nx)) 
    {
      damp->x[j] = 1.0f;
    }
    else if (j < p->nb) 
    {
      int d = p->nb - j;
      damp->x[j] = exp(-(p->factor * d) * (p->factor * d));
    }
    else 
    {
      int d = j - (p->nb + p->nx - 1);
      damp->x[j] = exp(-(p->factor * d) * (p->factor * d));
    }
  }
}

static void get_damp_z(damping_t *damp, const config_t *p)
{
  for (int i = 0; i < p->nzz; i++) 
  {
    // same as nb <= i <= nz + nb
    if ((unsigned)(i - p->nb) < (p->nz)) 
    {
      damp->z[i] = 1.0f;
    }
    else if (i < p->nb) 
    {
      int d = p->nb - i;
      damp->z[i] = exp(-(p->factor * d) * (p->factor * d));
    }
    else 
    {
      int d = i - (p->nb + p->nz - 1);
      damp->z[i] = exp(-(p->factor * d) * (p->factor * d));
    }
   }
}

damping_t static get_damp(const config_t *p)
{
  damping_t damp;

  damp.x = (float *)calloc(p->nxx, sizeof(float));
  damp.z = (float *)calloc(p->nzz, sizeof(float));

  if (!damp.x || !damp.z) 
  {
    perror("Could not allocate damping");
    exit(EXIT_FAILURE);
  }

  get_damp_z(&damp, p);
  get_damp_x(&damp, p);

  return damp;
}

static void allocate_fields(const config_t *p, fields_t *fld)
{
  size_t n = p->nxx * p->nzz;

  fld->txx    = (float *)calloc(n, sizeof(float));
  fld->tzz    = (float *)calloc(n, sizeof(float));
  fld->txz    = (float *)calloc(n, sizeof(float));
  fld->vx     = (float *)calloc(n, sizeof(float));
  fld->vz     = (float *)calloc(n, sizeof(float));
  fld->calc_p = (float *)calloc(n, sizeof(float));
}

// TODO
void register_seismogram() {}

void fd(const config_t *p, model_t *m, fields_t *fld)
{
  allocate_fields(p, fld);
  set_boundary(p, m);

  damping_t damp = get_damp(p);

  int nxx = p->nxx;
  int nzz = p->nzz;

  float inv_dx = 1.0f / p->dx;
  float inv_dz = 1.0f / p->dz;

  float dt = p->dt;

  #pragma omp parallel
  {
    for (size_t t = 0; t < p->nt; t++) 
    {
      #pragma omp single
      inject_source(p, fld, t);

      fd_pressure8E2T(fld, m, nxx, nzz, inv_dx, inv_dz, dt, &damp);
      fd_velocity8E2T(fld, m, nxx, nzz, inv_dx, inv_dz, dt, &damp);

      #pragma omp single
      if (p->snap_bool)
          get_snapshots(p, fld, t);
    } 
  } 

  free(fld->txx); free(fld->tzz); free(fld->txz);
  free(fld->vx);  free(fld->vz);  free(fld->calc_p);
  free(damp.x);   free(damp.z);
}

