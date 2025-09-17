#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "bin.h"
#include "wavelet.h"
#include "par.h"
#include "fd.h"

#define OUTPUT_PATH "data/output/"

struct timespec start, end;

int main(void)
{
  clock_gettime(CLOCK_MONOTONIC, &start);

  model_par_t model_par;
  model_par.rho_path = "data/input/model_rho_2d_1150x648.bin";
  model_par.vp_path  = "data/input/model_vp_2d_1150x648.bin";
  model_par.vs_path  = "data/input/model_vs_2d_1150x648.bin";

  model_par.nx     = 1150;
  model_par.nz     = 648;
  model_par.nb     = 100;
  model_par.factor = 0.015f;
  model_par.dx     = 5.0f;
  model_par.dz     = 5.0f;

  model_par.nxx = model_par.nx + 2 * model_par.nb;
  model_par.nzz = model_par.nz + 2 * model_par.nb;

  wavelet_t wavelet;
  wavelet.nt      = 5001;
  wavelet.dt      = 4.4e-4f;
  wavelet.fmax    = 30.0f;
  wavelet.wavelet = ricker(wavelet.nt, wavelet.dt, wavelet.fmax);

  snap_t snap;
  snap.snap_bool  = 1;
  snap.snap_num   = 100;
  snap.snap_ratio = wavelet.nt / snap.snap_num;

  model_t model;
  model.vp  = malloc(sizeof(float) * model_par.nx * model_par.nz);
  model.vs  = malloc(sizeof(float) * model_par.nx * model_par.nz);
  model.rho = malloc(sizeof(float) * model_par.nx * model_par.nz);

  read2D(model_par.vp_path,  model.vp, sizeof(float), model_par.nx, model_par.nz);
  read2D(model_par.vs_path,  model.vs, sizeof(float), model_par.nx, model_par.nz);
  read2D(model_par.rho_path, model.rho, sizeof(float), model_par.nx, model_par.nz);

  fields_t fld = {0};

  fd(&model_par, &wavelet, &snap, &model, &fld);

  clock_gettime(CLOCK_MONOTONIC, &end);

  double elapsed = (end.tv_sec - start.tv_sec)
                  + (end.tv_nsec - start.tv_nsec) / 1e9;

  printf("Elapsed: %.4f seconds\n", elapsed);

  free(m.vp);
  free(m.vs);
  free(m.rho);
  free(wavelet.wavelet);

  return 0;
}

