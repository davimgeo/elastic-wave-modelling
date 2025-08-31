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
      .rho_path = "data/input/homo_model/model_rho_2d_501x501.bin",
      .vp_path  = "data/input/homo_model/model_vp_2d_501x501.bin",
      .vs_path  = "data/input/homo_model/model_vs_2d_501x501.bin",
      .nx       = 501,
      .nz       = 501,
      .dx       = 10.0f,  
      .dz       = 10.0f
  };

  waveletPar wpar = 
  {
      .nt   = 1001,
      .dt   = 1e-3f,
      .fmax = 30
  };
  
  snapshots snap = 
  {
      .snap_bool = 1,
      .snap_num  = 5
  };

  geomPar geom = 
  {
      .sIdx = 250,
      .sIdz = 250
  };

  fdFields fld = {0};

  wpar.wavelet = ricker(wpar.nt, wpar.dt, wpar.fmax);

  mpar.vp  = read_f32_bin_model(mpar.vp_path, mpar.nx, mpar.nz);
  write_f32_bin_model("data/output/vp.bin", mpar.vp, mpar.nx, mpar.nz);
  mpar.vs  = read_f32_bin_model(mpar.vs_path, mpar.nx, mpar.nz);
  mpar.rho = read_f32_bin_model(mpar.rho_path, mpar.nx, mpar.nz);

  fd(&fld, mpar.vp, mpar.vs, mpar.rho, mpar.nx, mpar.nz, wpar.nt,
    wpar.wavelet, wpar.dt, mpar.dx, mpar.dz, geom.sIdx, geom.sIdz, &snap);

  char *filenames[] = {"txx.bin", "txz.bin", "tzz.bin", "vx.bin", "vz.bin"};
  char paths[5][256]; 

  for (int i = 0; i < 5; i++) 
  {
      snprintf(paths[i], sizeof(paths[i]), "%s%s", OUTPUT_PATH, filenames[i]);
  }

  // write_f32_bin_model(txx_path, fld.txx, mpar.nx, mpar.nz);
  // write_f32_bin_model(txz_path, fld.txz, mpar.nx, mpar.nz);
  // write_f32_bin_model(tzz_path, fld.tzz, mpar.nx, mpar.nz);
  // write_f32_bin_model(vx_path,  fld.vx,  mpar.nx, mpar.nz);
  // write_f32_bin_model(vz_path,  fld.vz,  mpar.nx, mpar.nz);

  clock_gettime(CLOCK_MONOTONIC, &end);

  double elapsed = (end.tv_sec - start.tv_sec)
                  + (end.tv_nsec - start.tv_nsec) / 1e9;

  printf("Elapsed: %.4f seconds\n", elapsed);
}

