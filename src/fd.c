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

  /* pad top and bottom respectively */
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

/* its truncating in 648~ for some reason */
static inline void generate_p(fdFields *fld, float *calc_p, int nx, int nz)
{
	for (size_t i = 0; i < nz; i++) 
	{
    for (size_t j = 0; j < nx; j++) 
    {
      calc_p[i + j * nz] = 0.5f * (fld->txx[i + j * nz] + fld->tzz[i + j * nz]);
    }
	}
}

static void 
get_snapshots
( int       snap_ratio, 
  fdFields *fld, 
  float    *calc_vp,
  int       time_step,
  int       nx, 
  int       nz )
{
  if (!(time_step % snap_ratio))
  {
    generate_p(fld, calc_vp, nx, nz);

    const char *filenames[] = {
      "vp_%dx%d_tid_%d.bin",
      "vx_%dx%d_tid_%d.bin", 
      "vz_%dx%d_tid_%d.bin"
    };

    float *fields[] = { calc_vp, fld->vx, fld->vz };

    char current_snap[256];
    char full_path[512];

    for (int i = 0; i < 3; i++) 
    {
      snprintf(current_snap, sizeof(current_snap), filenames[i], nx, nz, time_step);
      snprintf(full_path, sizeof(full_path), "%s%s", OUTPUT_PATH, current_snap);

      write_f32_bin_model(full_path, fields[i], nx, nz);
    }

    printf("Generating snapshot num %d...\n", time_step);
  }
}

static void 
fd_velocity_8E2T(modelPar *mdl, fdFields *fld, waveletPar *wav)
{
	float dx = mdl->dx;
	float dz = mdl->dz;
	float dt = wav->dt;
	int   nx = mdl->nx;
	int   nz = mdl->nz;

	float *rho = mdl->rho;

	#pragma omp parallel for schedule(static)
	for (int index = 0; index < nx * nz; index++)
	{
		int i = index % nz;
		int j = index / nz;

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

			float inv_rho_1 =
				1.0f / (0.5f * (rho[i + j * nz] + rho[i + (j + 1) * nz]));

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

			float inv_rho_2 =
				1.0f / (0.5f * (rho[i + j * nz] + rho[(i + 1) + j * nz]));

			fld->vz[i + j * nz] += dt * inv_rho_2 * (dtxz_dx + dtzz_dz);
		}
	}
}

static void 
fd_pressure_8E2T
( modelPar   *mdl, 
  fdFields   *fld, 
  waveletPar *wav )
{
	int   nx  = mdl->nx;
	int   nz  = mdl->nz;
	float dx  = mdl->dx;
	float dz  = mdl->dz;
	float dt  = wav->dt;

	float *vp  = mdl->vp;
	float *vs  = mdl->vs;
	float *rho = mdl->rho;

	#pragma omp parallel for schedule(static)
	for (int index = 0; index < nx * nz; index++)
	{
		int i = index % nz;
		int j = index / nz;

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

			float vp2     = vp[i + j * nz] * vp[i + j * nz];
			float vs2     = vs[i + j * nz] * vs[i + j * nz];

			float lambda  = rho[i + j * nz] * (vp2 - 2.0f * vs2);
			float mi      = rho[i + j * nz] * vs2;

			fld->txx[i + j * nz] += dt * ((lambda + 2.0f * mi) * dvx_dx +
								   lambda * dvz_dz);
			fld->tzz[i + j * nz] += dt * ((lambda + 2.0f * mi) * dvz_dz +
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

			float vs2       = vs[i       +       j * nz] * vs[      i + j * nz];
			float vs2_xp    = vs[(i + 1) +       j * nz] * vs[(i + 1) + j * nz];
			float vs2_zp    = vs[i       + (j + 1) * nz] * vs[i + (j + 1) * nz];
			float vs2_xp_zp = 
        vs[(i + 1) + (j + 1) * nz] * 
        vs[(i + 1) + (j + 1) * nz];

			float mi1 = rho[      i + j * nz]       * vs2;
			float mi2 = rho[(i + 1) + j * nz]       * vs2_xp;
			float mi3 = rho[i       + (j + 1) * nz] * vs2_zp;
			float mi4 = rho[(i + 1) + (j + 1) * nz] * vs2_xp_zp;

			float mi_avg = 4.0f / ((1.0f / mi1) + (1.0f / mi2) +
								   (1.0f / mi3) + (1.0f / mi4));

			fld->txz[i + j * nz] += dt * mi_avg * (dvx_dz + dvz_dx);
		}
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
  size_t n = mdl->nx * mdl->nz;

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

  float *calc_p = (float *)calloc(n, sizeof(float));

  if (!calc_p) 
  {
    perror("Could not allocate calculated vp\n");
    return;
  }

  int snap_ratio = wav->nt / snap->snap_num;

  // int nxx = mld->nx + 2 * mld->nb;
  // int nzz = mld->nz + 2 * mld->nb;

  for (int t = 0; t < wav->nt; t++)
  {
    int s_idx = geom->sIdz + geom->sIdx * mdl->nz;

    fld->txx[s_idx] += wav->wavelet[t] / (mdl->dx * mdl->dz);
    fld->tzz[s_idx] += wav->wavelet[t] / (mdl->dx * mdl->dz);

    fd_velocity_8E2T(mdl, fld, wav);
    fd_pressure_8E2T(mdl, fld, wav);

    // self.Vx *= self.damp2D
    // self.Vz *= self.damp2D
    // self.Txx *= self.damp2D
    // self.Tzz *= self.damp2D
    // self.Txz *= self.damp2D

    if (snap->snap_bool)
      get_snapshots(snap_ratio, fld, calc_p, t, mdl->nx, mdl->nz);
  }
}
