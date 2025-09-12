#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "bin.h"
#include "fd.h"
#include "par.h"

#define PI 3.1415926f

/* finite-difference coefficients */
#define FDM8E1 6.97545e-4f
#define FDM8E2 9.57031e-3f
#define FDM8E3 7.97526e-2f
#define FDM8E4 1.19628906f

#define OUTPUT_PATH "data/output/snapshots/"
     
/* declaring config struct globally inside this file */
static config_t *p = NULL;

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

void set_boundary()
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
      vp_ext[(i + nb) + (j + nb) * nzz]  = p->vp[i + j * nz];
      vs_ext[(i + nb) + (j + nb) * nzz]  = p->vs[i + j * nz];
      rho_ext[(i + nb) + (j + nb) * nzz] = p->rho[i + j * nz];
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
      vp_ext[i + j * nzz]  = vp_ext[i + nb * nzz];
      vs_ext[i + j * nzz]  = vs_ext[i + nb * nzz];
      rho_ext[i + j * nzz] = rho_ext[i + nb * nzz];

      vp_ext[i + (nx + nb + j) * nzz]  = vp_ext[i + (nx + nb - 1) * nzz];
      vs_ext[i + (nx + nb + j) * nzz]  = vs_ext[i + (nx + nb - 1) * nzz];
      rho_ext[i + (nx + nb + j) * nzz] = rho_ext[i + (nx + nb - 1) * nzz];  
    }
  }

  /* swap pointers to new arr */
  free(p->vp); p->vp = vp_ext;
  free(p->vs); p->vs = vs_ext;
  free(p->rho); p->rho = rho_ext;
}

void generate_p()
{
  int nxx = p->nxx;
  int nzz = p->nzz;

	for (size_t i = 0; i < nzz; i++) 
	{
    for (size_t j = 0; j < nxx; j++) 
    {
      p->calc_p[i + j * nzz] = 
        0.5f * (p->txx[i + j * nzz] + p->tzz[i + j * nzz]);
    }
	}
}

void get_snapshots(int time_step)
{
  if (!(time_step % p->snap_ratio))
  {
    generate_p();

    const char *filenames[] = {
      "p_%dx%d_tid_%d.bin",
      "vx_%dx%d_tid_%d.bin", 
      "vz_%dx%d_tid_%d.bin"
    };

    float *fields[] = { p->calc_p, p->vx, p->vz };

    char current_snap[256], full_path[512];

    #define FIELDS_SIZE 3
    for (int i = 0; i < FIELDS_SIZE; i++) 
    {
      snprintf(current_snap, sizeof(current_snap), filenames[i], p->nxx, p->nzz, time_step);
      snprintf(full_path, sizeof(full_path), "%s%s", OUTPUT_PATH, current_snap);

      write_f32_bin_model(full_path, fields[i], p->nxx, p->nzz);
    }

    printf("Generating snapshot num %d...\n", time_step);
  }
}

void apply_boundary(damping_t *damp)
{
  int nzz = p->nzz;
  int nxx = p->nxx;

  #pragma omp parallel for
  for (size_t index = 0; index < nxx * nzz; index++) 
  {
    int i = index % nzz;
    int j = index / nzz;

    p->calc_p[i + j * nzz] *= (damp->x[j] * damp->z[i]); 
    p->txz[i + j * nzz]    *= (damp->x[j] * damp->z[i]); 
    p->vx[i + j * nzz]     *= (damp->x[j] * damp->z[i]); 
    p->vz[i + j * nzz]     *= (damp->x[j] * damp->z[i]); 
  }
}

void fd_velocity_8E2T()
{
    const int nxx = p->nxx;
    const int nzz = p->nzz;

    const float inv_dx = 1.0f / p->dx;
    const float inv_dz = 1.0f / p->dz;
    const float dt = p->dt;

    float *restrict txx = p->txx;
    float *restrict txz = p->txz;
    float *restrict tzz = p->tzz;
    float *restrict vx  = p->vx;
    float *restrict vz  = p->vz;
    float *restrict rho = p->rho;

    #pragma omp parallel for
    for (int j = 0; j < nxx; ++j)
    {
      for (int i = 0; i < nzz; ++i)
      {
        if ((i >= 4) && (i < nzz - 4) && (j > 4) && (j < nxx - 4))
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

          float rho_inv = 1.0f / (0.5f * (rho[i + j * nzz] + rho[i + (j + 1) * nzz]));

          vx[i + j * nzz] += dt * rho_inv * (dtxx_dx + dtxz_dz);
        }

        if ((i > 4) && (i < nzz - 4) && (j >= 4) && (j < nxx - 4))
        {
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

            float rho_inv = 1.0f / (0.5f * (rho[i + j * nzz] + rho[(i + 1) + j * nzz]));

            vz[i + j * nzz] += dt * rho_inv * (dtxz_dx + dtzz_dz);
        }
    }
  }
}


void fd_pressure_8E2T()
{
    const int nxx = p->nxx;
    const int nzz = p->nzz;

    const float inv_dx = 1.0f / p->dx;
    const float inv_dz = 1.0f / p->dz;
    const float dt = p->dt;

    float *restrict vx  = p->vx;
    float *restrict vz  = p->vz;
    float *restrict txx = p->txx;
    float *restrict tzz = p->tzz;
    float *restrict txz = p->txz;
    float *restrict vp  = p->vp;
    float *restrict vs  = p->vs;
    float *restrict rho = p->rho;

    #pragma omp parallel for
    for (int j = 0; j < nxx; ++j)
    {
      for (int i = 0; i < nzz; ++i)
      {
        if ((i >= 4) && (i < nzz - 4) && (j >= 4) && (j < nxx - 4))
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

            float vp2 = vp[i + j * nzz] * vp[i + j * nzz];
            float vs2 = vs[i + j * nzz] * vs[i + j * nzz];

            float lamb = rho[i + j * nzz] * (vp2 - 2.0f * vs2);
            float mi   = rho[i + j * nzz] * vs2;

            txx[i + j * nzz] += dt * ((lamb + 2.0f * mi) * dvx_dx + lamb * dvz_dz);
            tzz[i + j * nzz] += dt * ((lamb + 2.0f * mi) * dvz_dz + lamb * dvx_dx);
        }

        if ((i >= 4) && (i <= nzz - 4) && (j >= 4) && (j <= nxx - 4))
        {
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

            float vs2       = vs[i + j * nzz] * vs[i + j * nzz];
            float vs2_xp    = vs[(i + 1) + j * nzz] * vs[(i + 1) + j * nzz];
            float vs2_zp    = vs[i + (j + 1) * nzz] * vs[i + (j + 1) * nzz];
            float vs2_xp_zp = vs[(i + 1) + (j + 1) * nzz] * vs[(i + 1) + (j + 1) * nzz];

            float mi1 = rho[i + j * nzz]         * vs2;
            float mi2 = rho[(i + 1) + j * nzz]   * vs2_xp;
            float mi3 = rho[i + (j + 1) * nzz]   * vs2_zp;
            float mi4 = rho[(i + 1) + (j + 1) * nzz] * vs2_xp_zp;

            float mi_avg = 4.0f / ((1.0f / mi1) + (1.0f / mi2) + (1.0f / mi3) + (1.0f / mi4));

            txz[i + j * nzz] += dt * mi_avg * (dvx_dz + dvz_dx);
      }
    }
  }
}


void inject_source(size_t t)
{
  int s_idx = (p->sIdz + p->nb) + (p->sIdx + p->nb) * p->nzz;

  p->txx[s_idx] += p->wavelet[t] / (p->dx * p->dz);
  p->tzz[s_idx] += p->wavelet[t] / (p->dx * p->dz);
}

damping_t get_damp()
{
  damping_t damp;

  damp.x = (float *)calloc(p->nxx, sizeof(float));
  damp.z = (float *)calloc(p->nzz, sizeof(float));

  if (!damp.x || !damp.z) 
  {
    perror("Could not allocate damping");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < p->nzz; i++) 
  {
    if ((i >= p->nb) && (i < p->nz + p->nb)) 
    {
      damp.z[i] = 1.0f;
    }
    else if (i < p->nb) 
    {
      int d = p->nb - i;
      damp.z[i] = exp(-(p->factor * d) * (p->factor * d));
    }
    else 
    {
      int d = i - (p->nb + p->nz - 1);
      damp.z[i] = exp(-(p->factor * d) * (p->factor * d));
    }
  }

  for (int j = 0; j < p->nxx; j++) 
  {
    if ((j >= p->nb) && (j < p->nx + p->nb)) 
    {
      damp.x[j] = 1.0f;
    }
    else if (j < p->nb) 
    {
      int d = p->nb - j;
      damp.x[j] = exp(-(p->factor * d) * (p->factor * d));
    }
    else 
    {
      int d = j - (p->nb + p->nx - 1);
      damp.x[j] = exp(-(p->factor * d) * (p->factor * d));
    }
  }

  return damp;
}

void allocate_fields()
{
  size_t n = p->nxx * p->nzz;

  p->txx = (float *)calloc(n, sizeof(float));
  p->tzz = (float *)calloc(n, sizeof(float));
  p->txz = (float *)calloc(n, sizeof(float));
  p->vx  = (float *)calloc(n, sizeof(float));
  p->vz  = (float *)calloc(n, sizeof(float));

  p->calc_p = (float *)calloc(n, sizeof(float));

  if (!p->txx || !p->tzz || !p->txz || 
      !p->vx || !p->vz || !p->calc_p)
  {
    perror("Could not allocate fields\n");
    return;
  }
}

void register_seismogram() {}

void fd(config_t *config)
{
  p = config;

  allocate_fields();

  set_boundary();

  write_f32_bin_model("data/output/vp.bin", p->vp, p->nxx, p->nzz);

  damping_t damp = get_damp();  

  for (size_t t = 0; t < p->nt; t++)
  {
    //register_seismogram();

    inject_source(t);

    fd_velocity_8E2T();
    fd_pressure_8E2T();

    apply_boundary(&damp);

    if (p->snap_bool)
      get_snapshots(t);
  }

  free(p->txx); free(p->tzz);
  free(p->txz); free(p->vx);
  free(p->calc_p); 
  free(damp.x); free(damp.z);
}
