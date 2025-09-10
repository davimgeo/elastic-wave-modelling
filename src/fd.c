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

void set_boundary(config_t *p)
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

void generate_p(config_t *p)
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

void get_snapshots(config_t *p, int time_step)
{
  if (!(time_step % p->snap_ratio))
  {
    generate_p(p);

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

void apply_boundary(config_t *p, damping_t *damp)
{
  int nzz = p->nzz;
  int nxx = p->nxx;

	#pragma omp parallel for schedule(static)
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

void fd_velocity_8E2T(config_t *p)
{
	int nxx = p->nxx;
	int nzz = p->nzz;

	#pragma omp parallel for schedule(static)
	for (int index = 0; index < nxx * nzz; index++)
	{
		int i = index % nzz;
		int j = index / nzz;

		if ((i >= 4) && (i < nzz - 4) && (j > 4) && (j < nxx - 4))
		{
			float dtxx_dx =
				(FDM8E1 * (p->txx[i + (j - 4) * nzz] - p->txx[i + (j + 3) * nzz]) +
				 FDM8E2 * (p->txx[i + (j + 2) * nzz] - p->txx[i + (j - 3) * nzz]) +
				 FDM8E3 * (p->txx[i + (j - 2) * nzz] - p->txx[i + (j + 1) * nzz]) +
				 FDM8E4 * (p->txx[i + j * nzz] - p->txx[i + (j - 1) * nzz])) / p->dx;

			float dtxz_dz =
				(FDM8E1 * (p->txz[(i - 3) + j * nzz] - p->txz[(i + 4) + j * nzz]) +
				 FDM8E2 * (p->txz[(i + 3) + j * nzz] - p->txz[(i - 2) + j * nzz]) +
				 FDM8E3 * (p->txz[(i - 1) + j * nzz] - p->txz[(i + 2) + j * nzz]) +
				 FDM8E4 * (p->txz[(i + 1) + j * nzz] - p->txz[i + j * nzz])) / p->dz;

			float rho_inv = 1.0f / (0.5f * (p->rho[i + j * nzz] + p->rho[i + (j + 1) * nzz]));

			p->vx[i + j * nzz] += p->dt * rho_inv * (dtxx_dx + dtxz_dz);
		}

		if ((i > 4) && (i < nzz - 4) && (j >= 4) && (j < nxx - 4))
		{
			float dtxz_dx =
				(FDM8E1 * (p->txz[i + (j - 3) * nzz] - p->txz[i + (j + 4) * nzz]) +
				 FDM8E2 * (p->txz[i + (j + 3) * nzz] - p->txz[i + (j - 2) * nzz]) +
				 FDM8E3 * (p->txz[i + (j - 1) * nzz] - p->txz[i + (j + 2) * nzz]) +
				 FDM8E4 * (p->txz[i + (j + 1) * nzz] - p->txz[i + j * nzz])) / p->dx;

			float dtzz_dz =
				(FDM8E1 * (p->tzz[(i - 4) + j * nzz] - p->tzz[(i + 3) + j * nzz]) +
				 FDM8E2 * (p->tzz[(i + 2) + j * nzz] - p->tzz[(i - 3) + j * nzz]) +
				 FDM8E3 * (p->tzz[(i - 2) + j * nzz] - p->tzz[(i + 1) + j * nzz]) +
				 FDM8E4 * (p->tzz[i + j * nzz] - p->tzz[(i - 1) + j * nzz])) / p->dz;

			float rho_inv = 1.0f / (0.5f * (p->rho[i + j * nzz] + p->rho[(i + 1) + j * nzz]));

			p->vz[i + j * nzz] += p->dt * rho_inv * (dtxz_dx + dtzz_dz);
		}
	}
}

void fd_pressure_8E2T(config_t *p)
{
	int   nxx  = p->nxx;
	int   nzz  = p->nzz;

	#pragma omp parallel for schedule(static)
	for (int index = 0; index < nxx * nzz; index++)
	{
		int i = index % nzz;
		int j = index / nzz;

		if ((i >= 4) && (i < nzz - 4) && (j >= 4) && (j < nxx - 4))
		{
			float dvx_dx =
				(FDM8E1 * (p->vx[i + (j - 3) * nzz] - p->vx[i + (j + 4) * nzz]) +
				 FDM8E2 * (p->vx[i + (j + 3) * nzz] - p->vx[i + (j - 2) * nzz]) +
				 FDM8E3 * (p->vx[i + (j - 1) * nzz] - p->vx[i + (j + 2) * nzz]) +
				 FDM8E4 * (p->vx[i + (j + 1) * nzz] - p->vx[i + j * nzz])) / p->dx;

			float dvz_dz =
				(FDM8E1 * (p->vz[(i - 3) + j * nzz] - p->vz[(i + 4) + j * nzz]) +
				 FDM8E2 * (p->vz[(i + 3) + j * nzz] - p->vz[(i - 2) + j * nzz]) +
				 FDM8E3 * (p->vz[(i - 1) + j * nzz] - p->vz[(i + 2) + j * nzz]) +
				 FDM8E4 * (p->vz[(i + 1) + j * nzz] - p->vz[i + j * nzz])) / p->dz;

			float vp2 = p->vp[i + j * nzz] * p->vp[i + j * nzz];
			float vs2 = p->vs[i + j * nzz] * p->vs[i + j * nzz];

			float lamb = p->rho[i + j * nzz] * (vp2 - 2.0f * vs2);
			float mi   = p->rho[i + j * nzz] * vs2;

			p->txx[i + j * nzz] += p->dt * ((lamb + 2.0f * mi) * dvx_dx + lamb * dvz_dz);
			p->tzz[i + j * nzz] += p->dt * ((lamb + 2.0f * mi) * dvz_dz + lamb * dvx_dx);
		}

		if ((i >= 4) && (i <= nzz - 4) && (j >= 4) && (j <= nxx - 4))
		{
			float dvx_dz =
				(FDM8E1 * (p->vx[(i - 4) + j * nzz] - p->vx[(i + 3) + j * nzz]) +
				 FDM8E2 * (p->vx[(i + 2) + j * nzz] - p->vx[(i - 3) + j * nzz]) +
				 FDM8E3 * (p->vx[(i - 2) + j * nzz] - p->vx[(i + 1) + j * nzz]) +
				 FDM8E4 * (p->vx[i + j * nzz] - p->vx[(i - 1) + j * nzz])) / p->dz;

			float dvz_dx =
				(FDM8E1 * (p->vz[i + (j - 4) * nzz] - p->vz[i + (j + 3) * nzz]) +
				 FDM8E2 * (p->vz[i + (j + 2) * nzz] - p->vz[i + (j - 3) * nzz]) +
				 FDM8E3 * (p->vz[i + (j - 2) * nzz] - p->vz[i + (j + 1) * nzz]) +
				 FDM8E4 * (p->vz[i + j * nzz] - p->vz[i + (j - 1) * nzz])) / p->dx;

			float vs2       = p->vs[i + j * nzz] * p->vs[i + j * nzz];
			float vs2_xp    = p->vs[(i + 1) + j * nzz] * p->vs[(i + 1) + j * nzz];
			float vs2_zp    = p->vs[i + (j + 1) * nzz] * p->vs[i + (j + 1) * nzz];
			float vs2_xp_zp = p->vs[(i + 1) + (j + 1) * nzz] * p->vs[(i + 1) + (j + 1) * nzz];

			float mi1 = p->rho[i + j * nzz] * vs2;
			float mi2 = p->rho[(i + 1) + j * nzz] * vs2_xp;
			float mi3 = p->rho[i + (j + 1) * nzz] * vs2_zp;
			float mi4 = p->rho[(i + 1) + (j + 1) * nzz] * vs2_xp_zp;

			float mi_avg = 4.0f / ((1.0f / mi1) + (1.0f / mi2) + (1.0f / mi3) + (1.0f / mi4));

			p->txz[i + j * nzz] += p->dt * mi_avg * (dvx_dz + dvz_dx);
		}
	}
}

void inject_source(config_t *p, size_t t)
{
  int s_idx = (p->sIdz + p->nb) + (p->sIdx + p->nb) * p->nzz;

  p->txx[s_idx] += p->wavelet[t] / (p->dx * p->dz);
  p->tzz[s_idx] += p->wavelet[t] / (p->dx * p->dz);
}

damping_t get_damp(config_t *p)
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

void allocate_fields(config_t *p)
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

void fd(config_t *p)
{
  allocate_fields(p);

  set_boundary(p);

  damping_t damp = get_damp(p);  

  for (size_t t = 0; t < p->nt; t++)
  {
    //register_seismogram();

    inject_source(p, t);

    fd_velocity_8E2T(p);
    fd_pressure_8E2T(p);

    apply_boundary(p, &damp);

    if (p->snap_bool)
      get_snapshots(p, t);
  }

  free(p->txx); free(p->tzz);
  free(p->txz); free(p->vx);
  free(p->calc_p); 
  free(damp.x); free(damp.z);
}
