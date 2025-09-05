#include <time.h>

#include "../debug/debug.h"

#include "par.h"
#include "bin.h"
#include "fd.h"

#define OUTPUT_PATH "data/output/"

int main(void)
{
  struct timespec start, end;
  clock_gettime(CLOCK_MONOTONIC, &start);

  modelPar mpar = 
  {
      .rho_path = "data/input/model_rho_2d_1150x648.bin",
      .vp_path  = "data/input/model_vp_2d_1150x648.bin",
      .vs_path  = "data/input/model_vs_2d_1150x648.bin",
      .nx       = 1150,
      .nz       = 648,
      .nb       = 30,
      .dx       = 5.0f,  
      .dz       = 5.0f
  };

  waveletPar wpar = 
  {
      .nt   = 5001,
      .dt   = 4.4e-4f,
      .fmax = 20
  };
  
  snapshots snap = 
  {
      .snap_bool = 1,
      .snap_num  = 30
  };

  geomPar geom = 
  {
      .sIdx = 575,
      .sIdz = 83
  };

  fdFields fld = {0};

  wpar.wavelet = ricker(wpar.nt, wpar.dt, wpar.fmax);

  mpar.vp  = read_f32_bin_model(mpar.vp_path, mpar.nx, mpar.nz);
  mpar.vs  = read_f32_bin_model(mpar.vs_path, mpar.nx, mpar.nz);
  mpar.rho = read_f32_bin_model(mpar.rho_path, mpar.nx, mpar.nz);

  //set_boundary(&fld, &mpar);

  //int nxx = mpar.nx + 2 * mpar.nb;
  //int nzz = mpar.nz + 2 * mpar.nb;

  //write_f32_bin_model("data/output/vp.bin", mpar.vp, nxx, nzz);
  //write_f32_bin_model("data/output/vs.bin", mpar.vs, nxx, nzz);
  //write_f32_bin_model("data/output/rho.bin", mpar.rho, nxx, nzz);

  fd(&fld, &mpar, &wpar, &geom, &snap);

  clock_gettime(CLOCK_MONOTONIC, &end);

  double elapsed = (end.tv_sec - start.tv_sec)
                  + (end.tv_nsec - start.tv_nsec) / 1e9;

  printf("Elapsed: %.4f seconds\n", elapsed);
}

