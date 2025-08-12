#include <stdlib.h>

#include "debug/debug.h"

#include "src/par.h"
#include "src/read_bin.h"
#include "src/wavelet.h"
#include "src/fd.h"

int main(int argc, char** argv)
{
    modelPar mpar = {
       .path = "data/input/model_vp_2d_1150x648.bin",
        .nx  = 1150,
        .nz  = 648
    };

    waveletPar wpar = {
        .nt   = 1001,
        .dt   = 0.004f,
        .fmax = 10
    };

    float *ricker_wavelet = ricker(wpar.nt, wpar.dt, wpar.fmax);

    // f32arr_to_txt("data/ricker_nt1001_dt4e-3.txt", ricker_wavelet, wpar.nt);

    float *model = read_f32_bin_model(mpar.path, mpar.nx, mpar.nz);

    // print_f32_2d_arr(model, msize.nx, msize.nz);

    free(model);
    return 0;
}

