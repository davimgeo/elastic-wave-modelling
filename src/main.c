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

  config_t c;

  c.rho_path = "data/input/model_rho_2d_1150x648.bin";
  c.vp_path  = "data/input/model_vp_2d_1150x648.bin";
  c.vs_path  = "data/input/model_vs_2d_1150x648.bin";

  c.nx     = 1150;
  c.nz     = 648;
  c.nb     = 100;
  c.factor = 0.015f;
  c.dx     = 5.0f;
  c.dz     = 5.0f;

  c.nxx = c.nx + 2 * c.nb;
  c.nzz = c.nz + 2 * c.nb;

  c.nt      = 5001;
  c.dt      = 4.4e-4f;
  c.fmax    = 30.0f;
  c.wavelet = ricker(c.nt, c.dt, c.fmax);

  c.snap_bool  = true;
  c.snap_num   = 100;
  c.snap_ratio = c.nt / c.snap_num;

  c.sIdx = 0.5 * (c.nxx - c.nb);
  c.sIdz = 0;

  c.vp  = read_f32_bin_model(c.vp_path,  c.nx, c.nz);
  c.vs  = read_f32_bin_model(c.vs_path,  c.nx, c.nz);
  c.rho = read_f32_bin_model(c.rho_path, c.nx, c.nz);

  fd(&c);

  clock_gettime(CLOCK_MONOTONIC, &end);

  double elapsed = (end.tv_sec - start.tv_sec)
                  + (end.tv_nsec - start.tv_nsec) / 1e9;

  printf("Elapsed: %.4f seconds\n", elapsed);
  return 0;
}

