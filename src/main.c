#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "io_bin.h"
#include "fd_wavelet.h"
#include "fd_core.h" 
#include "cfg_par.h"

#define OUTPUT_PATH "data/output/"

int main(void)
{
  struct timespec start, end;
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

  cfg.snap_bool  = 1;
  cfg.snap_num   = 100;
  cfg.snap_ratio = cfg.nt / cfg.snap_num;

  model_t model;
  model.vp  = malloc(sizeof(float) * cfg.nx * cfg.nz);
  model.vs  = malloc(sizeof(float) * cfg.nx * cfg.nz);
  model.rho = malloc(sizeof(float) * cfg.nx * cfg.nz);

  read2D(cfg.vp_path,  model.vp,  sizeof(float), cfg.nx, cfg.nz);
  read2D(cfg.vs_path,  model.vs,  sizeof(float), cfg.nx, cfg.nz);
  read2D(cfg.rho_path, model.rho, sizeof(float), cfg.nx, cfg.nz);

  cfg.src_f_lines = 2;
  cfg.src_f_cols  = 2;

  cfg.src_x = malloc(sizeof(float) * cfg.src_f_lines * cfg.src_f_cols);
  cfg.src_z = malloc(sizeof(float) * cfg.src_f_lines * cfg.src_f_cols);

  cfg.r_f_lines = 1151;
  cfg.r_f_cols  = 2;

  cfg.rcv_x = malloc(sizeof(float) * cfg.r_f_lines * cfg.r_f_cols);
  cfg.rcv_z = malloc(sizeof(float) * cfg.r_f_lines * cfg.r_f_cols);

  read_receivers("data/input/receivers.txt",
                  cfg.r_f_lines, cfg.r_f_cols,
                  cfg.rcv_x, cfg.rcv_z);

  read_sources("data/input/sources.txt",
                cfg.src_f_lines, cfg.src_f_cols,
                cfg.src_x, cfg.src_z);

  fields_t fld = {0};
  fd(&cfg, &fld, &model);

  clock_gettime(CLOCK_MONOTONIC, &end);
  double elapsed = (end.tv_sec - start.tv_sec)
                  + (end.tv_nsec - start.tv_nsec) / 1e9;

  printf("Elapsed: %.4f seconds\n", elapsed);

  free(model.vp);
  free(model.vs);
  free(model.rho);
  free(cfg.wavelet);
  free(cfg.src_x);
  free(cfg.src_z);
  free(cfg.rcv_x);
  free(cfg.rcv_z);

  return 0;
}

