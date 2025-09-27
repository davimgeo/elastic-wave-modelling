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

#define FREE(ptr) do { free(ptr); ptr = NULL; } while(0)

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

  for (int i = 0; i < p->nzz; i++) 
  {
    // same as nb <= i <= nz + nb
    if ((unsigned)(i - p->nb) < (p->nz)) 
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
    // same as nb <= j <= nz + nb
    if ((unsigned)(j - p->nb) < (p->nx)) 
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

void get_snapshots(const config_t *p, const fields_t *fld, int time_step)
{
  if (!(time_step % p->snap_ratio)) 
  {
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

void fd(const config_t *p, model_t *m, fields_t *fld)
{
  allocate_fields(p, fld);
  set_boundary(p, m);

  write2D("data/output/vp.bin", m->vp, sizeof(float), p->nxx, p->nzz);

  int seismogram_size = p->nt * p->r_f_lines;
  float* seismogram = (float *)calloc((size_t)seismogram_size, sizeof(float));
  if (!seismogram) { perror("Could not allocate Seismogram\n"); exit(EXIT_FAILURE); }

  damping_t damp = get_damp(p);

  int nxx = p->nxx;
  int nzz = p->nzz;
  float inv_dx = 1.0f / p->dx;
  float inv_dz = 1.0f / p->dz;
  float inv_dx_dz = inv_dx * inv_dz;
  float dt = p->dt;

  for (size_t s = 0; s < p->src_f_lines - 1; s++)
  {
    int sIdx = p->src_x[s];
    int sIdz = p->src_z[s];

    int s_idx = (sIdz + p->nb) + (sIdx + p->nb) * p->nzz;

    for (size_t t = 0; t < p->nt; t++)
    {
      /*
       * inject source
       */
      fld->txx[s_idx] += p->wavelet[t] * inv_dx_dz;
      fld->tzz[s_idx] += p->wavelet[t] * inv_dx_dz;

      /*
       * pressure update
       */
      #pragma omp parallel for schedule(static)
      for (int j = 4; j < nxx - 4; j++)
      {
        for (int i_ = 4; i_ < nzz - 4; i_++)
        {
          float dvx_dx =
            (FDM8E1 * (fld->vx[i_ + (j - 3) * nzz] - fld->vx[i_ + (j + 4) * nzz]) +
             FDM8E2 * (fld->vx[i_ + (j + 3) * nzz] - fld->vx[i_ + (j - 2) * nzz]) +
             FDM8E3 * (fld->vx[i_ + (j - 1) * nzz] - fld->vx[i_ + (j + 2) * nzz]) +
             FDM8E4 * (fld->vx[i_ + (j + 1) * nzz] - fld->vx[i_ + j * nzz])) * inv_dx;

          float dvz_dz =
            (FDM8E1 * (fld->vz[(i_ - 3) + j * nzz] - fld->vz[(i_ + 4) + j * nzz]) +
             FDM8E2 * (fld->vz[(i_ + 3) + j * nzz] - fld->vz[(i_ - 2) + j * nzz]) +
             FDM8E3 * (fld->vz[(i_ - 1) + j * nzz] - fld->vz[(i_ + 2) + j * nzz]) +
             FDM8E4 * (fld->vz[(i_ + 1) + j * nzz] - fld->vz[i_ + j * nzz])) * inv_dz;

          float dvx_dz =
            (FDM8E1 * (fld->vx[(i_ - 4) + j * nzz] - fld->vx[(i_ + 3) + j * nzz]) +
             FDM8E2 * (fld->vx[(i_ + 2) + j * nzz] - fld->vx[(i_ - 3) + j * nzz]) +
             FDM8E3 * (fld->vx[(i_ - 2) + j * nzz] - fld->vx[(i_ + 1) + j * nzz]) +
             FDM8E4 * (fld->vx[i_ + j * nzz] - fld->vx[(i_ - 1) + j * nzz])) * inv_dz;

          float dvz_dx =
            (FDM8E1 * (fld->vz[i_ + (j - 4) * nzz] - fld->vz[i_ + (j + 3) * nzz]) +
             FDM8E2 * (fld->vz[i_ + (j + 2) * nzz] - fld->vz[i_ + (j - 3) * nzz]) +
             FDM8E3 * (fld->vz[i_ + (j - 2) * nzz] - fld->vz[i_ + (j + 1) * nzz]) +
             FDM8E4 * (fld->vz[i_ + j * nzz] - fld->vz[i_ + (j - 1) * nzz])) * inv_dx;

          float vp2 = m->vp[i_ + j * nzz] * m->vp[i_ + j * nzz];
          float vs2 = m->vs[i_ + j * nzz] * m->vs[i_ + j * nzz];
          float vs2_xp = m->vs[(i_ + 1) + j * nzz] * m->vs[(i_ + 1) + j * nzz];
          float vs2_zp = m->vs[i_ + (j + 1) * nzz] * m->vs[i_ + (j + 1) * nzz];
          float vs2_xp_zp = m->vs[(i_ + 1) + (j + 1) * nzz] * m->vs[(i_ + 1) + (j + 1) * nzz];

          float lamb = m->rho[i_ + j * nzz] * (vp2 - 2.0f * vs2);
          float mi   = m->rho[i_ + j * nzz] * vs2;

          float mi1 = m->rho[i_ + j * nzz] * vs2;
          float mi2 = m->rho[(i_ + 1) + j * nzz] * vs2_xp;
          float mi3 = m->rho[i_ + (j + 1) * nzz] * vs2_zp;
          float mi4 = m->rho[(i_ + 1) + (j + 1) * nzz] * vs2_xp_zp;
          float mi_avg = 4.0f / ((1.0f / mi1) + (1.0f / mi2) + (1.0f / mi3) + (1.0f / mi4));

          fld->txx[i_ + j * nzz] += dt * ((lamb + 2.0f * mi) * dvx_dx + lamb * dvz_dz);
          fld->tzz[i_ + j * nzz] += dt * ((lamb + 2.0f * mi) * dvz_dz + lamb * dvx_dx);
          fld->txz[i_ + j * nzz] += dt * mi_avg * (dvx_dz + dvz_dx);

          float damp_prod = damp.x[j] * damp.z[i_];
          fld->txx[i_ + j * nzz] *= damp_prod;
          fld->tzz[i_ + j * nzz] *= damp_prod;
          fld->txz[i_ + j * nzz] *= damp_prod;

          fld->calc_p[i_ + j * nzz] = 0.5f * 
            (fld->txx[i_ + j * nzz] + fld->tzz[i_ + j * nzz]);
        }
      }

      /*
       * velocity update
       */
      #pragma omp parallel for schedule(static)
      for (int j = 4; j < nxx - 4; j++)
      {
        for (int i_ = 4; i_ < nzz - 4; i_++)
        {
          float dtxx_dx =
            (FDM8E1 * (fld->txx[i_ + (j - 4) * nzz] - fld->txx[i_ + (j + 3) * nzz]) +
             FDM8E2 * (fld->txx[i_ + (j + 2) * nzz] - fld->txx[i_ + (j - 3) * nzz]) +
             FDM8E3 * (fld->txx[i_ + (j - 2) * nzz] - fld->txx[i_ + (j + 1) * nzz]) +
             FDM8E4 * (fld->txx[i_ + j * nzz] - fld->txx[i_ + (j - 1) * nzz])) * inv_dx;

          float dtxz_dz =
            (FDM8E1 * (fld->txz[(i_ - 3) + j * nzz] - fld->txz[(i_ + 4) + j * nzz]) +
             FDM8E2 * (fld->txz[(i_ + 3) + j * nzz] - fld->txz[(i_ - 2) + j * nzz]) +
             FDM8E3 * (fld->txz[(i_ - 1) + j * nzz] - fld->txz[(i_ + 2) + j * nzz]) +
             FDM8E4 * (fld->txz[(i_ + 1) + j * nzz] - fld->txz[i_ + j * nzz])) * inv_dz;

          float dtxz_dx =
            (FDM8E1 * (fld->txz[i_ + (j - 3) * nzz] - fld->txz[i_ + (j + 4) * nzz]) +
             FDM8E2 * (fld->txz[i_ + (j + 3) * nzz] - fld->txz[i_ + (j - 2) * nzz]) +
             FDM8E3 * (fld->txz[i_ + (j - 1) * nzz] - fld->txz[i_ + (j + 2) * nzz]) +
             FDM8E4 * (fld->txz[i_ + (j + 1) * nzz] - fld->txz[i_ + j * nzz])) * inv_dx;

          float dtzz_dz =
            (FDM8E1 * (fld->tzz[(i_ - 4) + j * nzz] - fld->tzz[(i_ + 3) + j * nzz]) +
             FDM8E2 * (fld->tzz[(i_ + 2) + j * nzz] - fld->tzz[(i_ - 3) + j * nzz]) +
             FDM8E3 * (fld->tzz[(i_ - 2) + j * nzz] - fld->tzz[(i_ + 1) + j * nzz]) +
             FDM8E4 * (fld->tzz[i_ + j * nzz] - fld->tzz[(i_ - 1) + j * nzz])) * inv_dz;

          float rho_inv  = 1.0f / (0.5f * (m->rho[i_ + j * nzz] + m->rho[i_ + (j + 1) * nzz]));
          float rho_inv2 = 1.0f / (0.5f * (m->rho[i_ + j * nzz] + m->rho[(i_ + 1) + j * nzz]));

          fld->vx[i_ + j * nzz] += dt * rho_inv  * (dtxx_dx + dtxz_dz);
          fld->vz[i_ + j * nzz] += dt * rho_inv2 * (dtxz_dx + dtzz_dz);

          float damp_prod = damp.x[j] * damp.z[i_];
          fld->vx[i_ + j * nzz] *= damp_prod;
          fld->vz[i_ + j * nzz] *= damp_prod;
        }
      }

      if (p->snap_bool) { get_snapshots(p, fld, t); }
     
      /*
       * register seismogram
       */
      #pragma omp parallel for schedule(static)
      for (int i = 0; i < p->r_f_lines; i++)
      {
        int ix = (int)p->rcv_x[i] + p->nb;
        int iz = (int)p->rcv_z[i] + p->nb;

        int idx = iz + ix * p->nzz;

        seismogram[t + i * p->nt] = fld->calc_p[idx];
      }
    }
  }

  write2D("data/output/seismogram_txx_1150x648.bin", 
      seismogram, sizeof(float), p->nt, p->r_f_lines);

  FREE(seismogram);
  FREE(damp.x);
  FREE(damp.z);

  FREE(fld->txx);
  FREE(fld->tzz);
  FREE(fld->txz);
  FREE(fld->vx);
  FREE(fld->vz);
  FREE(fld->calc_p);
}



