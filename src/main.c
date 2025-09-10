#include <time.h>
#include <stdio.h>

#include "par.h"
#include "bin.h"
#include "fd.h"

#define OUTPUT_PATH "data/output/"

struct timespec start, end;

int main(void)
{
  clock_gettime(CLOCK_MONOTONIC, &start);

  config_t p;

  p.rho_path = "data/input/model_rho_2d_1150x648.bin";
  p.vp_path  = "data/input/model_vp_2d_1150x648.bin";
  p.vs_path  = "data/input/model_vs_2d_1150x648.bin";

  p.nx     = 1150;
  p.nz     = 648;
  p.nb     = 100;
  p.factor = 0.015f;
  p.dx     = 5.0f;
  p.dz     = 5.0f;

  p.nxx = p.nx + 2 * p.nb;
  p.nzz = p.nz + 2 * p.nb;

  p.nt      = 5001;
  p.dt      = 4.4e-4f;
  p.fmax    = 30.0f;
  p.wavelet = ricker(p.nt, p.dt, p.fmax);

  p.snap_bool  = true;
  p.snap_num   = 100;
  p.snap_ratio = p.nt / p.snap_num;

  p.sIdx = 0.5 * (p.nxx - p.nb);
  p.sIdz = 0;

  p.vp  = read_f32_bin_model(p.vp_path,  p.nx, p.nz);
  p.vs  = read_f32_bin_model(p.vs_path,  p.nx, p.nz);
  p.rho = read_f32_bin_model(p.rho_path, p.nx, p.nz);

  fd(&p);

  clock_gettime(CLOCK_MONOTONIC, &end);

  double elapsed = (end.tv_sec - start.tv_sec)
                  + (end.tv_nsec - start.tv_nsec) / 1e9;

  printf("Elapsed: %.4f seconds\n", elapsed);

  return 0;
}

