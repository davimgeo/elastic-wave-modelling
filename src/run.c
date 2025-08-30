#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../debug/debug.h"

#include "par.h"
#include "bin.h"
#include "fd.h"

#define OUTPUT_PATH "../data/output/"

int main (void)
{
  clock_t begin = clock();

  modelPar mpar = {
      .rho_path = "../data/input/model_rho_2d_1150x648.bin",
      .vp_path  = "../data/input/model_vp_2d_1150x648.bin",
      .vs_path  = "../data/input/model_vs_2d_1150x648.bin",
      .nx       = 501,
      .nz       = 501,
      .dx       = 10.0f,  
      .dz       = 10.0f
    };

  waveletPar wpar = {
      .nt   = 501,
      .dt   = 1e-3f,
      .fmax = 30
    };
  
  snapshots snap = {
      .snap_bool = 1,
      .snap_num  = 15
    };

  geomPar geom = {
      .sIdx = 250,
      .sIdz = 250
    };

  fdFields fld = {0};

  wpar.wavelet = ricker(wpar.nt, wpar.dt, wpar.fmax);

  fd(&fld, &mpar, &geom, &snap, &wpar);

  char txx_path[256], txz_path[256], tzz_path[256], vx_path[256], vz_path[256];

  snprintf(txx_path, sizeof(txx_path), "%s%s", OUTPUT_PATH, "txx.bin");
  snprintf(txz_path, sizeof(txz_path), "%s%s", OUTPUT_PATH, "txz.bin");
  snprintf(tzz_path, sizeof(tzz_path), "%s%s", OUTPUT_PATH, "tzz.bin");
  snprintf(vx_path,  sizeof(vx_path),  "%s%s", OUTPUT_PATH, "vx.bin");
  snprintf(vz_path,  sizeof(vz_path),  "%s%s", OUTPUT_PATH, "vz.bin");

  write_f32_bin_model(txx_path, fld.txx, mpar.nx, mpar.nz);
  // write_f32_bin_model(txz_path, fld.txz, mpar.nx, mpar.nz);
  // write_f32_bin_model(tzz_path, fld.tzz, mpar.nx, mpar.nz);
  // write_f32_bin_model(vx_path,  fld.vx,  mpar.nx, mpar.nz);
  // write_f32_bin_model(vz_path,  fld.vz,  mpar.nx, mpar.nz);

  clock_t end = clock();

  printf("Elapsed: %f seconds\n", (double)(end - begin) / (CLOCKS_PER_SEC * 1000.0f));
}

