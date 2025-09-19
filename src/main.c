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

  cfg.rho_path = "data/input/ro_351x1701_10m.bin";
  cfg.vp_path  = "data/input/vp_351x1701_10m.bin";
  cfg.vs_path  = "data/input/vs_351x1701_10m.bin";

  cfg.nx     = 1701;
  cfg.nz     = 351;
  cfg.nb     = 100;
  cfg.factor = 0.015f;
  cfg.dx     = 10.0f;
  cfg.dz     = 10.0f;

  cfg.nxx = cfg.nx + 2 * cfg.nb;
  cfg.nzz = cfg.nz + 2 * cfg.nb;

  cfg.nt      = 5001;
  cfg.dt      = 1e-3;
  cfg.fmax    = 15.0f;
  cfg.wavelet = ricker(cfg.nt, cfg.dt, cfg.fmax);

  cfg.snap_bool  = 1;
  cfg.snap_num   = 100;
  cfg.snap_ratio = cfg.nt / cfg.snap_num;

  model_t m;

  m.vp  = malloc(sizeof(float) * cfg.nx * cfg.nz);
  m.vs  = malloc(sizeof(float) * cfg.nx * cfg.nz);
  m.rho = malloc(sizeof(float) * cfg.nx * cfg.nz);

  read2D(cfg.vp_path, m.vp,  sizeof(float), cfg.nx, cfg.nz);
  read2D(cfg.vs_path, m.vs,  sizeof(float), cfg.nx, cfg.nz);
  read2D(cfg.rho_path, m.rho, sizeof(float), cfg.nx, cfg.nz);

  cfg.src_f_lines = 2;
  cfg.src_f_cols  = 2;

  cfg.src_x = malloc(sizeof(float) * cfg.src_f_lines * cfg.src_f_cols);
  cfg.src_z = malloc(sizeof(float) * cfg.src_f_lines * cfg.src_f_cols);

  cfg.r_f_lines = 172;
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

  fd(&cfg, &m, &fld);

  clock_gettime(CLOCK_MONOTONIC, &end);

  double elapsed = (end.tv_sec - start.tv_sec)
                  + (end.tv_nsec - start.tv_nsec) / 1e9;

  printf("Elapsed: %.4f seconds\n", elapsed);

  free(m.vp);
  free(m.vs);
  free(m.rho);
  free(cfg.wavelet);
  free(cfg.src_x);
  free(cfg.src_z);
  free(cfg.rcv_x);
  free(cfg.rcv_z);
  return 0;
}

