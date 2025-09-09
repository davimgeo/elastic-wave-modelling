#include <stdio.h>
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

void set_boundary(fdFields *fld, modelPar *mld)
{
  int nb  = mld->nb;
  int nx  = mld->nx;
  int nz  = mld->nz;

  int nxx = nx + 2 * nb;   
  int nzz = nz + 2 * nb;  

  size_t n = nxx * nzz;

  float *vp_ext  = (float *)calloc(n, sizeof(float));
  float *vs_ext  = (float *)calloc(n, sizeof(float));
  float *rho_ext = (float *)calloc(n, sizeof(float));

  /* copy original arr into ext */
  for (int j = 0; j < nx; j++) 
  {
    for (int i = 0; i < nz; i++) 
    {
      vp_ext[(i + nb) + (j + nb) * nzz]  = mld->vp[i + j * nz];
      vs_ext[(i + nb) + (j + nb) * nzz]  = mld->vs[i + j * nz];
      rho_ext[(i + nb) + (j + nb) * nzz] = mld->rho[i + j * nz];
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
  free(mld->vp); mld->vp = vp_ext;
  free(mld->vs); mld->vs = vs_ext;
  free(mld->rho); mld->rho = rho_ext;
}

void static inline generate_p(fdFields *fld, int nxx, int nzz)
{
	for (size_t i = 0; i < nzz; i++) 
	{
    for (size_t j = 0; j < nxx; j++) 
    {
      fld->calc_p[i + j * nzz] = 
        0.5f * (fld->txx[i + j * nzz] + fld->tzz[i + j * nzz]);
    }
	}
}

static void 
get_snapshots
( int       snap_ratio, 
  fdFields *fld, 
  int       time_step,
  int       nxx, 
  int       nzz )
{
  if (!(time_step % snap_ratio))
  {
    generate_p(fld, nxx, nzz);

    const char *filenames[] = {
      "p_%dx%d_tid_%d.bin",
      "vx_%dx%d_tid_%d.bin", 
      "vz_%dx%d_tid_%d.bin"
    };

    float *fields[] = { fld->calc_p, fld->vx, fld->vz };

    char current_snap[256];
    char full_path[512];

    const int FIELDS_SIZE = 3;
    for (int i = 0; i < FIELDS_SIZE; i++) 
    {
      snprintf(current_snap, sizeof(current_snap), filenames[i], nxx, nzz, time_step);
      snprintf(full_path, sizeof(full_path), "%s%s", OUTPUT_PATH, current_snap);

      write_f32_bin_model(full_path, fields[i], nxx, nzz);
    }

    printf("Generating snapshot num %d...\n", time_step);
  }
}

static void 
apply_boundary(fdFields *fld, float *damp1d, int nxx, int nzz)
{
  for (int i = 0; i < nzz; i++)
  {
      for (int j = 0; j < nxx; j++)
      {
        fld->vx[i + j * nzz] *= damp1d[i];
        fld->vz[i + j * nzz] *= damp1d[i];
        fld->txx[i + j * nzz] *= damp1d[i];
        fld->tzz[i + j * nzz] *= damp1d[i];
        fld->txz[i + j * nzz] *= damp1d[i];
        
        if (i < 1.0f) {
          printf("%f", damp1d[i]);
        }
      }
  }
}

static void 
fd_velocity_8E2T(modelPar *mdl, fdFields *fld, waveletPar *wav, float* damp1d)
{
	float dx = mdl->dx;
	float dz = mdl->dz;
	float dt = wav->dt;
	int   nxx = mdl->nx + 2 * mdl->nb;
	int   nzz = mdl->nz + 2 * mdl->nb;

	float *rho = mdl->rho;

	#pragma omp parallel for schedule(static)
	for (int index = 0; index < nxx * nzz; index++)
	{
		int i = index % nzz;
		int j = index / nzz;

		if ((i >= 4) && (i < nzz - 4) && (j > 4) && (j < nxx - 4))
		{
			float dtxx_dx =
				(FDM8E1 * (fld->txx[i + (j - 4) * nzz] - fld->txx[i + (j + 3) * nzz]) +
				 FDM8E2 * (fld->txx[i + (j + 2) * nzz] - fld->txx[i + (j - 3) * nzz]) +
				 FDM8E3 * (fld->txx[i + (j - 2) * nzz] - fld->txx[i + (j + 1) * nzz]) +
				 FDM8E4 * (fld->txx[i + j * nzz] - fld->txx[i + (j - 1) * nzz])) / dx;

			float dtxz_dz =
				(FDM8E1 * (fld->txz[(i - 3) + j * nzz] - fld->txz[(i + 4) + j * nzz]) +
				 FDM8E2 * (fld->txz[(i + 3) + j * nzz] - fld->txz[(i - 2) + j * nzz]) +
				 FDM8E3 * (fld->txz[(i - 1) + j * nzz] - fld->txz[(i + 2) + j * nzz]) +
				 FDM8E4 * (fld->txz[(i + 1) + j * nzz] - fld->txz[i + j * nzz])) / dz;

			float inv_rho_1 =
				1.0f / (0.5f * (rho[i + j * nzz] + rho[i + (j + 1) * nzz]));

			fld->vx[i + j * nzz] += dt * inv_rho_1 * (dtxx_dx + dtxz_dz);
		}

		if ((i > 4) && (i < nzz - 4) && (j >= 4) && (j < nxx - 4))
		{
			float dtxz_dx =
				(FDM8E1 * (fld->txz[i + (j - 3) * nzz] - fld->txz[i + (j + 4) * nzz]) +
				 FDM8E2 * (fld->txz[i + (j + 3) * nzz] - fld->txz[i + (j - 2) * nzz]) +
				 FDM8E3 * (fld->txz[i + (j - 1) * nzz] - fld->txz[i + (j + 2) * nzz]) +
				 FDM8E4 * (fld->txz[i + (j + 1) * nzz] - fld->txz[i + j * nzz])) / dx;

			float dtzz_dz =
				(FDM8E1 * (fld->tzz[(i - 4) + j * nzz] - fld->tzz[(i + 3) + j * nzz]) +
				 FDM8E2 * (fld->tzz[(i + 2) + j * nzz] - fld->tzz[(i - 3) + j * nzz]) +
				 FDM8E3 * (fld->tzz[(i - 2) + j * nzz] - fld->tzz[(i + 1) + j * nzz]) +
				 FDM8E4 * (fld->tzz[i + j * nzz] - fld->tzz[(i - 1) + j * nzz])) / dz;

			float inv_rho_2 =
				1.0f / (0.5f * (rho[i + j * nzz] + rho[(i + 1) + j * nzz]));

			fld->vz[i + j * nzz] += dt * inv_rho_2 * (dtxz_dx + dtzz_dz);
		}
	}
}

static void 
fd_pressure_8E2T(modelPar *mdl, fdFields *fld, waveletPar *wav, float *damp1d)
{
	int   nxx  = mdl->nx + 2 * mdl->nb;
	int   nzz  = mdl->nz + 2 * mdl->nb;
	float dx  = mdl->dx;
	float dz  = mdl->dz;
	float dt  = wav->dt;

	float *vp  = mdl->vp;
	float *vs  = mdl->vs;
	float *rho = mdl->rho;

	#pragma omp parallel for schedule(static)
	for (int index = 0; index < nxx * nzz; index++)
	{
		int i = index % nzz;
		int j = index / nzz;

		if ((i >= 4) && (i < nzz - 4) && (j >= 4) && (j < nxx - 4))
		{
			float dvx_dx =
				(FDM8E1 * (fld->vx[i + (j - 3) * nzz] - fld->vx[i + (j + 4) * nzz]) +
				 FDM8E2 * (fld->vx[i + (j + 3) * nzz] - fld->vx[i + (j - 2) * nzz]) +
				 FDM8E3 * (fld->vx[i + (j - 1) * nzz] - fld->vx[i + (j + 2) * nzz]) +
				 FDM8E4 * (fld->vx[i + (j + 1) * nzz] - fld->vx[i + j * nzz])) / dx;

			float dvz_dz =
				(FDM8E1 * (fld->vz[(i - 3) + j * nzz] - fld->vz[(i + 4) + j * nzz]) +
				 FDM8E2 * (fld->vz[(i + 3) + j * nzz] - fld->vz[(i - 2) + j * nzz]) +
				 FDM8E3 * (fld->vz[(i - 1) + j * nzz] - fld->vz[(i + 2) + j * nzz]) +
				 FDM8E4 * (fld->vz[(i + 1) + j * nzz] - fld->vz[i + j * nzz])) / dz;

			float vp2     = vp[i + j * nzz] * vp[i + j * nzz];
			float vs2     = vs[i + j * nzz] * vs[i + j * nzz];

			float lambda  = rho[i + j * nzz] * (vp2 - 2.0f * vs2);
			float mi      = rho[i + j * nzz] * vs2;

			fld->txx[i + j * nzz] += dt * ((lambda + 2.0f * mi) * dvx_dx +
								   lambda * dvz_dz);
			fld->tzz[i + j * nzz] += dt * ((lambda + 2.0f * mi) * dvz_dz +
								   lambda * dvx_dx);
		}

		if ((i >= 4) && (i <= nzz - 4) && (j >= 4) && (j <= nxx - 4))
		{
			float dvx_dz =
				(FDM8E1 * (fld->vx[(i - 4) + j * nzz] - fld->vx[(i + 3) + j * nzz]) +
				 FDM8E2 * (fld->vx[(i + 2) + j * nzz] - fld->vx[(i - 3) + j * nzz]) +
				 FDM8E3 * (fld->vx[(i - 2) + j * nzz] - fld->vx[(i + 1) + j * nzz]) +
				 FDM8E4 * (fld->vx[i + j * nzz] - fld->vx[(i - 1) + j * nzz])) / dz;

			float dvz_dx =
				(FDM8E1 * (fld->vz[i + (j - 4) * nzz] - fld->vz[i + (j + 3) * nzz]) +
				 FDM8E2 * (fld->vz[i + (j + 2) * nzz] - fld->vz[i + (j - 3) * nzz]) +
				 FDM8E3 * (fld->vz[i + (j - 2) * nzz] - fld->vz[i + (j + 1) * nzz]) +
				 FDM8E4 * (fld->vz[i + j * nzz] - fld->vz[i + (j - 1) * nzz])) / dx;

			float vs2       = vs[i       +       j * nzz] * vs[      i + j * nzz];
			float vs2_xp    = vs[(i + 1) +       j * nzz] * vs[(i + 1) + j * nzz];
			float vs2_zp    = vs[i       + (j + 1) * nzz] * vs[i + (j + 1) * nzz];
			float vs2_xp_zp = 
        vs[(i + 1) + (j + 1) * nzz] * 
        vs[(i + 1) + (j + 1) * nzz];

			float mi1 = rho[      i + j * nzz]       * vs2;
			float mi2 = rho[(i + 1) + j * nzz]       * vs2_xp;
			float mi3 = rho[i       + (j + 1) * nzz] * vs2_zp;
			float mi4 = rho[(i + 1) + (j + 1) * nzz] * vs2_xp_zp;

			float mi_avg = 4.0f / ((1.0f / mi1) + (1.0f / mi2) +
								   (1.0f / mi3) + (1.0f / mi4));

			fld->txz[i + j * nzz] += dt * mi_avg * (dvx_dz + dvz_dx);
		}
	}
}

void static 
inject_source
( fdFields   *fld,
  modelPar   *mdl, 
  geomPar    *geom, 
  waveletPar *wav, 
  size_t      t )
{
  int s_idx = geom->sIdz + geom->sIdx * mdl->nz;

  fld->txx[s_idx] += wav->wavelet[t] / (mdl->dx * mdl->dz);
  fld->tzz[s_idx] += wav->wavelet[t] / (mdl->dx * mdl->dz);
}

float* get_damp1D(modelPar *mdl)
{
  int nzz = mdl->nz + 2 * mdl->nb;

  float *damp1d = (float *)calloc(nzz, sizeof(float));

  if (!damp1d) 
  {
    perror("Could not allocate damp1d");
    return NULL;
  }

  for (int i = 0; i < nzz; i++) 
  {
    if ((i < mdl->nb)) 
    {
      int d = mdl->nb - i;
      damp1d[i] = exp(-(mdl->factor * d) * (mdl->factor * d));
    }
    else if ((i >= mdl->nb) && (i < mdl->nz)) 
    {
      damp1d[i] = 1.0f;
    }
    else 
    {
      int d = i - mdl->nz + 1;
      damp1d[i] = exp(-(mdl->factor * d) * (mdl->factor * d));
    }
    //printf("%f ", damp1d[i]);
  }

  return damp1d;
}

void static allocate_fields(fdFields *fld, int nxx, int nzz)
{
  size_t n = nxx * nzz;

  fld->txx = (float *)calloc(n, sizeof(float));
  fld->tzz = (float *)calloc(n, sizeof(float));
  fld->txz = (float *)calloc(n, sizeof(float));
  fld->vx  = (float *)calloc(n, sizeof(float));
  fld->vz  = (float *)calloc(n, sizeof(float));

  fld->calc_p = (float *)calloc(n, sizeof(float));

  if (!fld->txx || !fld->tzz || !fld->txz || 
      !fld->vx || !fld->vz || !fld->calc_p)
  {
    perror("Could not allocate fields\n");
    return;
  }
}

void 
fd
( fdFields   *fld, 
  modelPar   *mdl,
  waveletPar *wav,
  geomPar    *geom,
  snapshots  *snap )
{
  int nxx = mdl->nx + 2 * mdl->nb;
  int nzz = mdl->nz + 2 * mdl->nb;

  allocate_fields(fld, nxx, nzz);

  float *damp1d = get_damp1D(mdl);

  int snap_ratio = wav->nt / snap->snap_num;

  for (size_t t = 0; t < wav->nt; t++)
  {
    inject_source(fld, mdl, geom, wav, t);

    fd_velocity_8E2T(mdl, fld, wav, damp1d);
    fd_pressure_8E2T(mdl, fld, wav, damp1d);

    apply_boundary(fld, damp1d, nxx, nzz);

    if (snap->snap_bool)
      get_snapshots(snap_ratio, fld, t, nxx, nzz);
  }
}
