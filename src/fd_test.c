
#include <stdlib.h>

#include "fd.h"
# include <omp.h>
#include "bin.h"
#include "par.h"

#define SAFE_INV(x, eps) (1.0f / ((x) + (eps)))

#define OUTPUT_PATH "data/output/snapshots/"

#define PRINT(x) printf ("%f\n", (x))

#define IDX(i, j) ((i) + (j) * nz)

#define AVG(a, b) (0.5f * ((a) + (b)))

static void 
get_snapshots
( int    snap_ratio, 
  float *pressure, 
  int    time_step, 
  int    nx, 
  int    nz ) 
{
  if (!(time_step % snap_ratio)) 
      write_f32_bin_model(OUTPUT_PATH, pressure, nx, nz);
}

static void 
fd_velocity_2E2T
( modelPar *model, 
  fdFields *fld, 
  int       nx, 
  int       nz, 
  float     dx, 
  float     dz, 
  float     dt )
{
  #pragma omp parallel for
  for (int index = 0; index < nx*nz; index++)
  {
      int i = (int)(index % nz);
      int j = (int)(index / nz);

      if((i >= 0) && (i < nz-1) && (j > 0) && (j < nx)) 
      {
          
          float dtxx_dx = (Txx[i + j*nz]     - Txx[i + (j-1)*nz]) / dh;

          float dtxz_dz = (Txz[(i+1) + j*nz] - Txz[i + j*nz]) / dh;

          Vx[index] += dt*B*(dTxx_dx + dTxz_dz); 
      }

      if((i > 0) && (i < nz) && (j >= 0) && (j < nx-1)) 
      {
          
          float dTxz_dx = (txz[i + (j+1)*nz] - txz[i + j*nz]) / dh;

          float dTzz_dz = (tzz[i + j*nz]     - tzz[(i-1) + j*nz]) / dh;

          Vz[index] += dt*B*(dTxz_dx + dTzz_dz); 
      }            
}

static void 
fd_pressure_2E2T
( modelPar *model, 
  fdFields *fld, 
  int       nx, 
  int       nz, 
  float     dx, 
  float     dz, 
  float     dt )
{
  #pragma omp parallel for
  for (int index = 0; index < nx*nz; index++)
  {
      int i = (int)(index % nz);
      int j = (int)(index / nz);

      if((i >= 0) && (i < nz-1) && (j >= 0) && (j < nx-1)) 
      {    
          
          float dVx_dx = (FDM2E1*(Vx[i + (j+1)*nz] - Vx[i + j*nz])) / dh;

          float dVz_dz = (FDM2E1*(Vz[(i+1) + j*nz] - Vz[i + j*nz])) / dh;
                          
          Txx[index] += dt*((L + 2*M)*dVx_dx + L*dVz_dz);
          Tzz[index] += dt*((L + 2*M)*dVz_dz + L*dVx_dx);                    
      }

      if((i > 0) && (i < nz) && (j > 0) && (j < nx)) 
      {
          float dVx_dz = (FDM2E1*(Vx[i + j*nz]     - Vx[(i-1) + j*nz])) / dh;

          float dVz_dx = (FDM2E1*(Vz[i + j*nz]     - Vz[i + (j-1)*nz])) / dh;

          Txz[index] += dt*M*(dVx_dz + dVz_dx);
      }
  }
}

float* 
fd
( fdFields   *fld, 
  modelPar   *model, 
  geomPar    *geom, 
  snapshots  *snap, 
  waveletPar *wav )
{
  size_t n = model->nx * model->nz;

  /* initialize staggered grid arrays */
  fld->txx = (float *)calloc (n, sizeof (float));
  fld->tzz = (float *)calloc (n, sizeof (float));
  fld->txz = (float *)calloc (n, sizeof (float));
  fld->vx  = (float *)calloc (n, sizeof (float));
  fld->vz  = (float *)calloc (n, sizeof (float));

  if (!fld->txx || !fld->tzz || !fld->txz || !fld->vx || !fld->vz) 
      return NULL;
    
  int   nx = model->nx;
  int   nz = model->nz;
  float dx = model->dx;
  float dz = model->dz;

  int snap_ratio = wav->nt / snap->snap_num;

  for (int t = 0; t < wav->nt; t++) 
  {
      int s_idx = geom->sIdz + geom->sIdx * nz;

      fld->txx[s_idx] += wav->wavelet[t] / (dx * dz);
      fld->tzz[s_idx] += wav->wavelet[t] / (dx * dz);

      fd_velocity_2E2T (model, fld, nx, nz, dx, dz, wav->dt);
      fd_pressure_2E2T (model, fld, nx, nz, dx, dz, wav->dt);

      if (snap->snap_bool) 
          get_snapshots (snap_ratio, fld->txx, t, nx, nz);
        
    }

  free(fld->tzz);
  free(fld->txz);
  free(fld->vx);
  free(fld->vz);

  return fld->txx; 
}

