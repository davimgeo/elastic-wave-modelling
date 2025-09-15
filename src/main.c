#include <time.h>

#include "par.h"
#include "bin.h"
#include "wavelet.h"
#include "fd.h"

#define OUTPUT_PATH "data/output/"

struct timespec start, end;

int main(void)
{
  clock_gettime(CLOCK_MONOTONIC, &start);

  config_t cfg;

  cfg.rho_path = "data/input/model_rho_2d_1150x648.bin";
  cfg.vp_path  = "data/input/model_vp_2d_1150x648.bin";
  cfg.vs_path  = "data/input/model_vs_2d_1150x648.bin";

  cfg.nx     = 1150;
  cfg.nz     = 648;
  cfg.nb     = 100;
  cfg.factor = 0.015f;
  cfg.dx     = 5.0f;
  cfg.dz     = 5.0f;

  cfg.nxx = cfg.nx + 2 * cfg.nb;
  cfg.nzz = cfg.nz + 2 * cfg.nb;

  cfg.nt      = 5001;
  cfg.dt      = 4.4e-4f;
  cfg.fmax    = 30.0f;
  cfg.wavelet = ricker(cfg.nt, cfg.dt, cfg.fmax);

  cfg.snap_bool  = true;
  cfg.snap_num   = 100;
  cfg.snap_ratio = cfg.nt / cfg.snap_num;

  cfg.sIdx = 0.5 * (cfg.nxx - cfg.nb);
  cfg.sIdz = 0;

  model_t m;

  m.vp  = read_f32_bin_model(cfg.vp_path,  cfg.nx, cfg.nz);
  m.vs  = read_f32_bin_model(cfg.vs_path,  cfg.nx, cfg.nz);
  m.rho = read_f32_bin_model(cfg.rho_path, cfg.nx, cfg.nz);

  fields_t fld = {0};

  fd(&cfg, &m, &fld);

  clock_gettime(CLOCK_MONOTONIC, &end);

  double elapsed = (end.tv_sec - start.tv_sec)
                  + (end.tv_nsec - start.tv_nsec) / 1e9;

  printf("Elapsed: %.4f seconds\n", elapsed);
  return 0;
}

