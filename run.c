#include <stdlib.h>

#include "debug/debug.h"

#include "src/par.h"
#include "src/bin.h"
#include "src/wavelet.h"
#include "src/fd.h"

int main(void)
{
    modelPar mpar = {
        .rho_path = "data/input/model_rho_2d_1150x648.bin",
        .vp_path  = "data/input/model_vp_2d_1150x648.bin",
        .vs_path  = "data/input/model_vs_2d_1150x648.bin",
        .nx  = 1150,
        .nz  = 648,
        .dx  = 10.0f,  
        .dz  = 10.0f
    };

    waveletPar wpar = {
        .nt   = 2001,
        .dt   = 1e-3f,
        .fmax = 10
    };
    
    snapshots snap = {
        .snap_bool = 0,
        .snap_num  = 15,
    };

    geomPar geom = {
        .sIdx = 50,
        .sIdz = 550,
    };

    fdFields fld = {0};

    wpar.wavelet = ricker(wpar.nt, wpar.dt, wpar.fmax);

    mpar.rho = read_f32_bin_model(mpar.rho_path, mpar.nx, mpar.nz);
    mpar.vp  = read_f32_bin_model(mpar.vp_path,  mpar.nx, mpar.nz);
    mpar.vs  = read_f32_bin_model(mpar.vs_path,  mpar.nx, mpar.nz);
    //print_f32_2d_arr(mpar.vs, mpar.nx, mpar.nz);
    
    if (!mpar.rho || !mpar.vp || !mpar.vs) {
        free(mpar.rho); free(mpar.vp); free(mpar.vs);
        free(wpar.wavelet);
        return 1;
    }

    fld.txx = fd(&fld, &mpar, &geom, &snap, &wpar);
    print_f32_2d_arr(fld.txx, mpar.nx, mpar.nz);

    free(mpar.rho);
    free(mpar.vp);
    free(mpar.vs); 
    free(wpar.wavelet);

    free(fld.txx);
}
